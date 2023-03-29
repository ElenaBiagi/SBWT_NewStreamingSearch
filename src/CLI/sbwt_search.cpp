#include <string>
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SeqIO.hh"
#include "buffered_streams.hh"
#include "NewStreamingSearchRMQ.hh"
#include "NewStreamingSearchMin.hh"
#include "variants.hh"
#include "commands.hh"
#include <filesystem>
#include "suffix_group_optimization.hh"

using namespace std;

using namespace sbwt;



// Assumes values of v are -1 or larger
template <typename writer_t>
inline void print_vector(const vector<int64_t>& v, writer_t& out){
    // Fast manual integer-to-string conversion
    char buffer[32];
    char newline = '\n';
    for(int64_t x : v){
        int64_t i = 0;
        if(x == -1){
            buffer[0] = '1';
            buffer[1] = '-';
            i = 2;
        } else{
            while(x > 0){
                buffer[i++] = '0' + (x % 10);
                x /= 10;
            }
        }
        std::reverse(buffer, buffer + i);
        buffer[i] = ' ';
        out.write(buffer, i+1);
    }
    out.write(&newline, 1);
}

void print_LCS(const sdsl::int_vector<>& v,const string& fname) {
    std::ofstream csv_file(fname);
    // Write the contents of the int_vector to the CSV file
    for (size_t i = 0; i < v.size(); i++) {
        csv_file << v[i];
        if (i < v.size() - 1) {
            csv_file << ",";
        }
    }
    csv_file.close();
}


template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        int64_t t0 = cur_time_micros();
        vector<int64_t> out_buffer = sbwt.streaming_search(reader.read_buf, len);
        total_micros += cur_time_micros() - t0;

        number_of_queries += out_buffer.size();

        // Write out
        print_vector(out_buffer, writer);
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t new_run_queries_streaming_rmq(reader_t& reader, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const int64_t bit_len, const vector<int64_t>& C, const int64_t k, const sdsl::int_vector<>& LCS, const sdsl::rmq_succinct_sct<>& rmqLCS ){
    write_log("Inside NEW run streaming queries", LogLevel::MAJOR);
    int64_t new_total_micros = 0;
    int64_t new_number_of_queries = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        vector<pair<int64_t, int64_t>> new_search = new_streaming_search_rmq( DNA_bitvectors, DNA_rs, bit_len, C, k,LCS, rmqLCS,reader.read_buf, len);
        new_total_micros += cur_time_micros() - t0;
        new_number_of_queries += new_search.size();

        //print the vector
//        for (auto p :new_search){ cout << p.first << " , "<< p.second<< endl;}

        // test
//        vector<int64_t> out_buffer = sbwt.streaming_search(reader.read_buf, len);
//        const int64_t k = sbwt.get_k();
//        for (uint64_t q = 0; q < out_buffer.size();q++){
//            if (out_buffer[q]!=-1){
//                assert(new_search[q].first == k & new_search[q].second == out_buffer[q]);
//            }
//        }
    }
    write_log("us/query: " + to_string((double)new_total_micros / new_number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return new_number_of_queries;
}
template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t new_run_queries_streaming_min(reader_t& reader, const sbwt_t& sbwt, const sdsl::int_vector<>& LCS){

    int64_t new_total_micros = 0;
    int64_t new_number_of_queries = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = cur_time_micros();
        vector<pair<int64_t, int64_t>> new_search = new_streaming_search_min(sbwt, LCS, reader.read_buf, len);
        new_total_micros += cur_time_micros() - t0;
        new_number_of_queries += new_search.size();

        // test
//        vector<int64_t> out_buffer = sbwt.streaming_search(reader.read_buf, len);
//        const int64_t k = sbwt.get_k();
//        for (uint64_t q = 0; q < out_buffer.size();q++){
//            if (out_buffer[q]!=-1){
//                assert(new_search[q].first == k & new_search[q].second == out_buffer[q]);
//            }
//        }
    }
    write_log("us/query: " + to_string((double)new_total_micros / new_number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return new_number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_not_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t k = sbwt.get_k();
    vector<int64_t> out_buffer;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        for(int64_t i = 0; i < len - k + 1; i++){
            int64_t t0 = cur_time_micros();
            int64_t ans = sbwt.search(reader.read_buf + i);
            total_micros += cur_time_micros() - t0;
            number_of_queries++;
            out_buffer.push_back(ans);
        }

        print_vector(out_buffer, writer);
        out_buffer.clear();
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_file(const string& infile, const string& outfile, const sbwt_t& sbwt){
    reader_t reader(infile);
    writer_t writer(outfile);
    if(sbwt.has_streaming_query_support()){
        write_log("Running streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t new_run_file_rmq(const string& infile, const string& outfile, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const int64_t bit_len, const vector<int64_t>& C, const int64_t k, const sdsl::int_vector<>& LCS, const sdsl::rmq_succinct_sct<>& rmqLCS){
    reader_t reader(infile);
    writer_t writer(outfile);
    if(sbwt.has_streaming_query_support()){
        write_log("Running NEW streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return new_run_queries_streaming_rmq<sbwt_t, reader_t, writer_t>(reader, sbwt,DNA_bitvectors, DNA_rs, bit_len, C, k, LCS, rmqLCS);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t new_run_file_min(const string& infile, const string& outfile, const sbwt_t& sbwt, const sdsl::int_vector<>& LCS){
    reader_t reader(infile);
    writer_t writer(outfile);
    if(sbwt.has_streaming_query_support()){
        write_log("Running NEW streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return new_run_queries_streaming_min<sbwt_t, reader_t, writer_t>(reader, sbwt, LCS);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

// Returns number of queries executed
template<typename sbwt_t>
int64_t run_queries(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, bool gzip_output){

    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt);
        }
    }
    return n_queries_run;

}
// Returns number of queries executed
template<typename sbwt_t>
int64_t new_run_queries_rmq(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const int64_t bit_len, const vector<int64_t>& C, const int64_t k, const sdsl::int_vector<>& LCS, const sdsl::rmq_succinct_sct<>& rmqLCS, bool gzip_output){
    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += new_run_file_rmq<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, bit_len, C, k, LCS, rmqLCS);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += new_run_file_rmq<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, bit_len, C, k, LCS, rmqLCS);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += new_run_file_rmq<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, DNA_bitvectors, DNA_rs, bit_len, C, k, LCS, rmqLCS);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += new_run_file_rmq<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt,DNA_bitvectors, DNA_rs, bit_len, C, k, LCS, rmqLCS);
        }
    }
    return n_queries_run;

}

template<typename sbwt_t>
int64_t new_run_queries_min(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, const sdsl::int_vector<>& LCS, bool gzip_output){

    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in_gzip;
    typedef SeqIO::Reader<Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef Buffered_ofstream<zstr::ofstream> out_gzip;
    typedef Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = SeqIO::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += new_run_file_min<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, LCS);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += new_run_file_min<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, LCS);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += new_run_file_min<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, LCS);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += new_run_file_min<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, LCS);
        }
    }
    return n_queries_run;

}

std::vector<std::string> get_available_types(){
    return {"old", "new-rmq", "new-min"};
}

//lcs
void save_v(const std::string& filename, const sdsl::int_vector<>& v) {
    std::ofstream out(filename, std::ios::binary);
    sdsl::serialize(v, out);
    out.close();
}

void load_v(const std::string& filename, sdsl::int_vector<>& v) {
    std::ifstream in(filename, std::ios::binary);
    sdsl::load(v, in);
    in.close();
}
//rmq
void save_r(const std::string& filename, const sdsl::rmq_succinct_sct<>& v) {
    std::ofstream out(filename, std::ios::binary);
    sdsl::serialize(v, out);
    out.close();
}
void load_r(const std::string& filename, sdsl::rmq_succinct_sct<>& v){
    std::ifstream in(filename, std::ios::binary);
    sdsl::load(v, in);
    in.close();
}

int search_main(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all k-mers of all input reads.");

    vector<string> types = get_available_types();
    string all_types_string;
    for(string type : types) all_types_string += " " + type;

    options.add_options()
            ("o,out-file", "Output filename.", cxxopts::value<string>())
            ("i,index-file", "Index input file.", cxxopts::value<string>())
            ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
            ("z,gzip-output", "Writes output in gzipped form. This can shrink the output files by an order of magnitude.", cxxopts::value<bool>()->default_value("false"))
            ("t,type", "Decide which streaming search type you prefer. Available types: " + all_types_string,cxxopts::value<string>()->default_value("old")  )
            ("lcs", "Provide in input the LCS file if available.", cxxopts::value<string>()->default_value(""))
            ("rmq", "For the option --lcs new-rmq provide in input the rmqLCS file if available.", cxxopts::value<string>()->default_value(""))
            ("h,help", "Print usage")
            ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string type = opts["type"].as<string>();
    if(std::find(types.begin(), types.end(), type) == types.end()){
        cerr << "Error: unknown type: " << type << endl;
        cerr << "Available types are:" << all_types_string << endl;
        return 1;
    }

    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> input_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if(multi_file){
        input_files = readlines(queryfile);
    } else{
        input_files = {queryfile};
    }
    for(string file : input_files) check_readable(file);

    // Interpret output file
    string outfile = opts["out-file"].as<string>();
    bool gzip_output = opts["gzip-output"].as<bool>();
    vector<string> output_files;
    if(multi_file){
        output_files = readlines(outfile);
    } else{
        output_files = {outfile};
    }
    for(string file : output_files) check_writable(file);

    vector<string> variants = get_available_variants();

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    int64_t number_of_queries = 0;
    int64_t new_number_of_queries = 0;



    if (variant == "plain-matrix"){

        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        if (type == "old"){
            number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);

            int64_t total_micros = cur_time_micros() - micros_start;
            write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
        }
        else {
            const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
            const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
            const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
            const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;

            const sdsl::rank_support_v5<>& A_bits_rs = sbwt.get_subset_rank_structure().A_bits_rs;
            const sdsl::rank_support_v5<>& C_bits_rs = sbwt.get_subset_rank_structure().C_bits_rs;
            const sdsl::rank_support_v5<>& G_bits_rs = sbwt.get_subset_rank_structure().G_bits_rs;
            const sdsl::rank_support_v5<>& T_bits_rs = sbwt.get_subset_rank_structure().T_bits_rs;

            const int64_t k = sbwt.get_k();
            const vector<int64_t>& C = sbwt.get_C_array();

            vector<pair<int64_t,int64_t> > ans;
            const sdsl::bit_vector* DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
            const sdsl::rank_support_v5<>* DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};

            uint64_t  bit_len = A_bits.size(); // all the same

            string LCS_file = opts["lcs"].as<string>();
            if (LCS_file.empty()) {
                LCS_file = outfile + "LCS.sdsl";
                const sdsl::int_vector<> LCS = get_kmer_lcs(A_bits, C_bits, G_bits, T_bits, k);
                save_v(LCS_file, LCS);
            }
            sdsl::int_vector<> LCS;
            load_v(LCS_file,LCS);
            //print_LCS(LCS,outfile);

            if (type == "new-rmq"){
                string rmqLCS_file = opts["rmq"].as<string>();
                if (rmqLCS_file.empty()){
                    rmqLCS_file = outfile + "rmqLCS.sdsl";
                    const sdsl::rmq_succinct_sct<> rmqLCS(&LCS);
                    save_r(rmqLCS_file, rmqLCS);
                }
                sdsl::rmq_succinct_sct<> rmqLCS;
                load_r(rmqLCS_file, rmqLCS);
                new_number_of_queries += new_run_queries_rmq(input_files, output_files, sbwt, DNA_bitvectors, DNA_rs, bit_len, C, k, LCS, rmqLCS, gzip_output);
            }
            else if (type == "new-min"){
                new_number_of_queries += new_run_queries_min(input_files, output_files, sbwt, LCS, gzip_output);
            }
            int64_t new_total_micros = cur_time_micros() - micros_start;
            write_log("us/query end-to-end: " + to_string((double)new_total_micros / new_number_of_queries), LogLevel::MAJOR);
        }
    }
    else if (variant == "rrr-matrix"){
        rrr_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "mef-matrix"){
        mef_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "plain-split"){
        plain_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "rrr-split"){
        rrr_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "mef-split"){
        mef_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "plain-concat"){
        plain_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "mef-concat"){
        mef_concat_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "plain-subsetwt"){
        plain_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    else if (variant == "rrr-subsetwt"){
        rrr_sswt_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, gzip_output);
    }
    //int64_t total_micros = cur_time_micros() - micros_start;
    //write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);
    return 0;

}
