
std::vector<PileupLine> readFile(std::istream& in, bool parse_base_qualities, bool parse_mapping_qualities) {
    std::vector<PileupLine> result;
    if (in != cin) {
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
        result.reserve(lines.size());

        #pragma omp parallel for default(none) shared(s, lines, result)
        for (int i = 0; i < lines.size(); ++i) {
            PileupLine plp = parsePileupLine(lines[i])
            results.push_back(std::move(plp));
        }
        delete[] s;
    } else {
        // TODO read lines from stdin with posix getline()
    }
    return results;
}

std::vector<VCFRecord> callLikelihoodRatio(std::istream& in) {

}

std::vector<VCFRecord> callBayes(std::istream& in) {

}

std::vector<VCFRecord> callSiteMLError(std::istream& in, bool use_prior) {

}

std::vector<VCFRecord> callQualityBasedSimple(std::istream& in) {

}

std::vector<VCFRecord> callQualityBasedSamtools(std::istream& in) {

}
