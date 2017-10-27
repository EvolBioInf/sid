#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>

#include "call.hpp"
#include "lynch.hpp"
#include "stats.hpp"

std::vector<PileupLine> readFile(std::istream& in, const bool parse_base_qualities, const bool parse_mapping_qualities) {
	std::vector<PileupLine> result;
	for (std::string line; std::getline(in, line); ) {
		PileupLine plp = parsePileupLine(&line[0u], parse_base_qualities, parse_mapping_qualities);
		result.push_back(std::move(plp));
	}
	return result;
}

std::vector<PileupLine> readFileParallel(std::istream& in, const bool parse_base_qualities, const bool parse_mapping_qualities) {
	// read whole input file into buffer
	const auto start_pos = in.tellg();
	in.ignore(std::numeric_limits<std::streamsize>::max());
	const auto char_count = in.gcount();
	in.seekg(start_pos);
	char* s = new char[char_count+1];
	in.read(&s[0], char_count);
	s[char_count] = '\0';

	std::vector<char*> lines;
	char* saveptr;
	char* line = strtok_r(s, "\n", &saveptr);
	lines.push_back(line);
	while((line = strtok_r(nullptr, "\n", &saveptr)) != nullptr) {
		lines.push_back(line);
	}

	std::vector<PileupLine> result (lines.size());

	#pragma omp parallel for default(none) shared(s, lines, result)
	for (size_t i = 0; i < lines.size(); ++i) {
		PileupLine plp = parsePileupLine(lines[i], parse_base_qualities, parse_mapping_qualities);
		// result.push_back(std::move(plp));
		result[i] = std::move(plp);
	}
	delete[] s;
	return result;
}

std::vector<OutputRecord> callLikelihoodRatio(std::istream& in) {
	auto inputRecords = readFile(in, false, false);
	auto unique_profiles = countUniqueProfiles(inputRecords);

	unique_profiles.erase(
		std::remove_if(unique_profiles.begin(),
					   unique_profiles.end(),
					   [](const UniqueProfile& p) { return p.coverage < 4; }),
		unique_profiles.end()); 

	std::cerr << "# unique profiles: " << unique_profiles.size() << std::endl;

	// estimate genotype likelihoods
	auto nucleotide_distribution = computeNucleotideDistribution(unique_profiles);
	auto estimate = estimateProfileGenotypeLikelihoods(unique_profiles, nucleotide_distribution);

	std::cerr << std::scientific;
	std::cerr << "# heterozygosity: " << estimate.heterozygosity << std::endl;
	std::cerr << "# error: " << estimate.error_rate << std::endl;

	std::map<profile_t, size_t> index_of_profile;
	size_t i = 0;
	for (const auto& p : unique_profiles) {
		index_of_profile.emplace(p.profile, i++);
	}

	std::vector<Classification> profile_classes;
	profile_classes.reserve(unique_profiles.size());

	for (size_t i = 0; i < estimate.profile_likelihoods.size(); ++i) {
		double p1 = likelihoodRatioTest(estimate.profile_likelihoods[i].L_heterozygous,
										estimate.profile_likelihoods[i].L_homozygous);
		double p2 = likelihoodRatioTest(estimate.profile_likelihoods[i].L_homozygous,
										estimate.profile_likelihoods[i].L_heterozygous);

		std::string label {"hom"};
		if (p2 < 0.05) {
			label = "het";
		}
		profile_classes.push_back({
			label, "NN", static_cast<double>(p1), static_cast<double>(p2), "p_value", {}
		});
	}

	std::vector<OutputRecord> output_records;
	output_records.reserve(inputRecords.size());
	for (const auto& record : inputRecords) {
		auto index = index_of_profile.find(record.base_counts);
		if (index != index_of_profile.end()) {
			output_records.push_back({
				record.chromosome_name,
				record.position,
				profile_classes[index->second]
			});
		}
	}

	return output_records;
}

std::vector<OutputRecord> callBayes(std::istream& in) {
	auto inputRecords = readFile(in, false, false);
	auto unique_profiles = countUniqueProfiles(inputRecords);

	unique_profiles.erase(
		std::remove_if(unique_profiles.begin(),
					   unique_profiles.end(),
					   [](const UniqueProfile& p) { return p.coverage < 4; }),
		unique_profiles.end()); 

	std::cerr << "# unique profiles: " << unique_profiles.size() << std::endl;

	// estimate genotype likelihoods
	auto nucleotide_distribution = computeNucleotideDistribution(unique_profiles);
	auto estimate = estimateProfileGenotypeLikelihoods(unique_profiles, nucleotide_distribution);

	std::cerr << std::scientific;
	std::cerr << "# heterozygosity: " << estimate.heterozygosity << std::endl;
	std::cerr << "# error: " << estimate.error_rate << std::endl;

	// build reverse index
	std::map<profile_t, size_t> index_of_profile;
	size_t i = 0;
	for (const auto& p : unique_profiles) {
		index_of_profile.emplace(p.profile, i++);
	}

	// calculate genotypes
	std::vector<Classification> profile_classes;
	profile_classes.reserve(unique_profiles.size());

	for (size_t i = 0; i < estimate.profile_likelihoods.size(); ++i) {
		auto aposteriori_homozygous = estimate.profile_likelihoods[i].L_homozygous * (1 - estimate.heterozygosity);
		auto aposteriori_heterozygous = estimate.profile_likelihoods[i].L_heterozygous * estimate.heterozygosity;

		auto probability_homozygous = aposteriori_homozygous / (aposteriori_homozygous + aposteriori_heterozygous);
		auto probability_heterozygous = aposteriori_heterozygous / (aposteriori_homozygous + aposteriori_heterozygous);

		std::string label {"hom"};
		std::string genotype {"NN"};
		if (probability_heterozygous > probability_homozygous) {
			label = "het";
		}
		profile_classes.push_back({
			label, genotype, static_cast<double>(probability_homozygous), static_cast<double>(probability_heterozygous), "probability", {}
		});
	}

	// build output lines
	std::vector<OutputRecord> output_records;
	output_records.reserve(inputRecords.size());
	for (const auto& record : inputRecords) {
		auto index = index_of_profile.find(record.base_counts);
		if (index != index_of_profile.end()) {
			output_records.push_back({
				record.chromosome_name,
				record.position,
				profile_classes[index->second]
			});
		}
	}

	return output_records;
}

std::vector<OutputRecord> callSiteMLError(std::istream& in, const bool estimate_prior, double snp_prior, double error_threshold) {
	auto inputRecords = readFile(in, false, false);
	auto unique_profiles = countUniqueProfiles(inputRecords);

	std::map<profile_t, size_t> index_of_profile;
	size_t i = 0;
	for (const auto& p : unique_profiles) {
		index_of_profile.emplace(p.profile, i++);
	}

	if (estimate_prior) {
		auto unique_profiles_min_coverage = unique_profiles;
		unique_profiles_min_coverage.erase(
			std::remove_if(unique_profiles_min_coverage.begin(),
						   unique_profiles_min_coverage.end(),
						   [](const UniqueProfile& p) { return p.coverage < 4; }),
			unique_profiles_min_coverage.end()); 

		auto nucleotide_distribution = computeNucleotideDistribution(unique_profiles_min_coverage);
		auto estimate = estimateProfileGenotypeLikelihoods(unique_profiles_min_coverage, nucleotide_distribution);
		snp_prior = estimate.heterozygosity;
	}

	std::vector<Classification> profile_classes;
	profile_classes.reserve(unique_profiles.size());
	for (const UniqueProfile& p : unique_profiles) {
		std::array<int, 4> indices {0,1,2,3};
		std::sort(indices.begin(), indices.end(),
			[&p](const size_t a, const size_t b) {
				return p.profile[a] < p.profile[b];
			});
		// index of base with hightest (second highest) number of occurences on current site
		auto ref0 = indices[3];
		auto ref1 = indices[2];

		// homozygous: ML site error = (n2+n3+n4)/(n1+n2+n3+n4)
		double error1 = static_cast<double>(p.coverage - p.profile[ref0]) / static_cast<double>(p.coverage);
		if (error1 > error_threshold) {
			error1 = error_threshold;
		}
		long double l1 = homozygousLikelihood(p, error1, ref0);

		// heterozygous: ML site error = 1.5 * (n3+n4)/(n1+n2+n3+n4)
		double error2 = 1.5 * static_cast<double>(p.coverage - p.profile[ref0] - p.profile[ref1]) / static_cast<double>(p.coverage);
		if (error2 > error_threshold) {
			error2 = error_threshold;
		}
		long double l2 = heterozygousLikelihood(p, error2, ref0, ref1);

		if (snp_prior > 0) {
			l1 *= (1 - snp_prior);
			l2 *= snp_prior;
		}

		double p1 = likelihoodRatioTest(l2, l1);
		double p2 = likelihoodRatioTest(l1, l2);

		std::string label {"hom"};
		std::string genotype {"ACGT"[ref0], "ACGT"[ref0]};
		if (l2 > l1 && p2 < 0.05) {
			label = "het";
			genotype[1] = "ACGT"[ref1]; 
		}
		profile_classes.push_back({
			label, genotype, static_cast<double>(p1), static_cast<double>(p2), "p_value", {}
		});
	}
	std::vector<OutputRecord> output_records;
	output_records.reserve(inputRecords.size());
	for (const auto& record : inputRecords) {
		auto index = index_of_profile.find(record.base_counts);
		if (index != index_of_profile.end()) {
			output_records.push_back({
				record.chromosome_name,
				record.position,
				profile_classes[index->second]
			});
		}
	}

	return output_records;

}

std::vector<OutputRecord> callQualityBasedSimple(std::istream& in, const bool estimate_prior, double snp_prior) {
	auto input_records = readFile(in, true, true);
	std::vector<OutputRecord> output_records (input_records.size());

	if (estimate_prior) {
		auto unique_profiles = countUniqueProfiles(input_records);
		unique_profiles.erase(
			std::remove_if(unique_profiles.begin(),
						   unique_profiles.end(),
						   [](const UniqueProfile& p) { return p.coverage < 4; }),
			unique_profiles.end()); 

		auto nucleotide_distribution = computeNucleotideDistribution(unique_profiles);
		auto estimate = estimateProfileGenotypeLikelihoods(unique_profiles, nucleotide_distribution);
		snp_prior = estimate.heterozygosity;
	}


	#pragma omp parallel for default(none) shared(input_records, output_records, snp_prior)
	for (size_t i = 0; i < input_records.size(); ++i) {
		const PileupLine& plp = input_records[i];

		std::array<int, 4> indices {0,1,2,3};
		std::sort(indices.begin(), indices.end(),
			[&plp](const int a, const int b) {
				return plp.base_counts[a] < plp.base_counts[b];
			});
		auto ref0 = indices[3];
		auto ref1 = indices[2];

		// std::vector<double> errors;
		// errors.reserve(plp.bases.size());
		// std::transform(plp.base_qualities.begin(), plp.base_qualities.end(), plp.mapping_qualities.begin(),
		// 	std::back_inserter(errors), [](const uint8_t bq, const uint8_t mp) {
		// 		// convert from phred scale to probability
		// 		return pow(10., std::min(bq, mq) / -10.);
		// 	});
		long double log_probability_homozygous = 0;
		long double log_probability_heterozygous = 0;
		for (size_t j = 0; j < plp.bases.size(); ++j) {
			double error = pow(10., std::min(plp.base_qualities[j], plp.mapping_qualities[j]) / -10.);
			if (plp.bases[j] == "ACGT"[ref0]) {
				log_probability_homozygous += log(1 - error);
			} else {
				log_probability_homozygous += log(error);
			}
			if (plp.bases[j] == "ACGT"[ref0] || plp.bases[j] == "ACGT"[ref1]) {
				log_probability_heterozygous += log(1 - 2./3. * error);
			} else {
				log_probability_heterozygous += log(2./3. * error);
			}
		}

		auto logbinom = [](int n, int k) -> double {
			return log_gamma(n + 1) - log_gamma(n - k + 1) - log_gamma(k + 1);
		};
		int n = plp.base_counts[ref0] + plp.base_counts[ref1];
		int k = plp.base_counts[ref1];
		log_probability_heterozygous += logbinom(n, k) - n * logl(2); // ~ log({n choose k} * 1/2^n)

		// posterior probabilites
		long double pp1 = exp(log_probability_homozygous);
		long double pp2 = exp(log_probability_heterozygous);
		if (snp_prior > 0) {
			pp1 *= (1 - snp_prior);
			pp2 *= snp_prior;
		}

		double p1 = likelihoodRatioTest(pp2, pp1);
		double p2 = likelihoodRatioTest(pp1, pp2);

		std::string label {"hom"};
		std::string genotype {"ACGT"[ref0], "ACGT"[ref0]};
		if (p2 < 0.05) {
			label = "het";
			genotype[1] = "ACGT"[ref1]; 
		}
		output_records[i] = {plp.chromosome_name, plp.position,
			{ label, genotype, static_cast<double>(p1), static_cast<double>(p2), "p_value", {}}};
	}
	return output_records;
}
