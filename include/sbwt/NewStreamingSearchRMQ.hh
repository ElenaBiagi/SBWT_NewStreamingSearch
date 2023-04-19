#pragma once

#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include "SBWT.hh"
#include "variants.hh"
#include "sdsl/vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rmq_support.hpp"

#include "suffix_group_optimization.hh"

using namespace std;
using namespace sbwt;

// Already defined in SBWT.hh
int64_t get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Returns pairs (a_1, b_1), (a_2, b_2)..., such that
// - a_i is the length of the longest match ending at query[i]
// - b_i is the colex rank of one arbitrary k-mer having the longest match.
vector<pair<int64_t,int64_t> > new_streaming_search_rmq(const sdsl::bit_vector** DNA_bitvectors, const sdsl::rank_support_v5<>** DNA_rs, const vector<int64_t>& C, const int64_t k, const sdsl::int_vector<>& LCS, const sdsl::rmq_succinct_sct<>& rmqLCS, const char* input, const int64_t len){ //const plain_matrix_sbwt_t& sbwt
    // Search the first char
    // Start from the middle of the string
    vector<pair<int64_t,int64_t> > ans;

    const int64_t bit_len = DNA_bitvectors[0]->size();
    const int64_t start = bit_len/2;
    int64_t curr_pos = start;
    uint64_t match_len = 0;
    for(int64_t i = 0; i < len; i++){

        char c = static_cast<char>(input[i] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
        if(char_idx == -1) [[unlikely]] {// Invalid character
            match_len = 0;
            curr_pos = -1;
        } else {
            const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
            const sdsl::rank_support_v5<>& Bit_rs = *(DNA_rs[char_idx]);
            if (curr_pos == -1) [[unlikely]]{
                curr_pos = start; // restart from the middle of the string
            }
            // CASE 1
            if (Bit_v[curr_pos]) {
                match_len += (match_len < k);
                curr_pos = C[char_idx] + Bit_rs(curr_pos);
            } else {
                // CASE 2
                int64_t r = -1;
                int64_t r_index = curr_pos - 1;
                while (r_index > 0){
                    int word_len_r = std::min((int64_t)64, r_index);
                    int word_start_r = std::max((int64_t)0, (r_index-63));
                    uint64_t word_r = Bit_v.get_int(word_start_r, word_len_r);
                    if (word_r == 0) {
                        r_index -=word_len_r;
                    } else {
                        int leading_zeros= __builtin_clzll(word_r) - (64 - word_len_r); // count number of leading zeros
                        r = (r_index  - leading_zeros);
                        r = (r < 64) ? r -1 : r;
                        break;
                    }
                }

                int64_t l = bit_len;
                int64_t l_index = curr_pos + 1;
                while (l_index < bit_len){
                    int word_len_l = std::min(bit_len-l_index-2, (int64_t)64);
                    uint64_t word_l = Bit_v.get_int(l_index, word_len_l);
                    std::bitset<64> x(word_l);
                    if (word_l == 0) {
                        l_index +=word_len_l;
                    } else {
                        int trailing_zeros= __builtin_ctzll(word_l);
                        l = (l_index  + trailing_zeros);
                        break;
                    }
                }

                if (r >= 0) [[likely]]{
                    if (l < bit_len ) [[likely]] { // both ok
                        uint64_t r_lcs = LCS[rmqLCS(r + 1, curr_pos)];
                        uint64_t l_lcs = LCS[rmqLCS(curr_pos + 1, l)];

                        uint64_t max_lcs = max(r_lcs, l_lcs);
                        match_len = min(match_len, max_lcs) + 1;
                        int64_t next_dir = (max_lcs == r_lcs) ? r : l;
                        curr_pos = C[char_idx] + Bit_rs(next_dir);

                    } else { // r ok
                        uint64_t r_lcs = LCS[rmqLCS(r + 1, curr_pos)];
                        match_len = min(match_len, r_lcs) + 1;//((r_lcs <= match_len) ? r_lcs : match_len) + 1; // match_len -= -(match_len > l_lcs) & (match_len - l_lcs); // add ++
                        curr_pos = C[char_idx] + Bit_rs(r);
                    }
                } else if (l < bit_len) { // l ok
                    uint64_t l_lcs = LCS[rmqLCS(curr_pos + 1, l)];
                    match_len = min(match_len, l_lcs) + 1;//((l_lcs <= match_len) ? l_lcs : match_len)+1; // match_len -= -(match_len > l_lcs) & (match_len - l_lcs); // add ++
                    curr_pos = C[char_idx] + Bit_rs(l);
                } else [[unlikely]]{ // if the char is a DNA char but is not in the input SBWT
                    match_len = 0;
                    curr_pos = -1;
                }
            }
        }
        ans.push_back({match_len,curr_pos});
    }
    return ans;
}
