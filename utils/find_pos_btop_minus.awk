# Find reference position for a given query index based on the BTOP output from BLAST.
# If the query index refers to an insertion, i.e., there is no corresponding reference position,
# a fractional value is returned (e.g. 1000.5 if the insertion lies between reference nucleotides 1000 and 1001).
#
# Note: This script only works for BLAST hits on the MINUS STRAND strand of the reference genome!
#       The required input arguments differ from find_pos_btop_plus.awk (ref_start vs. ref_end).
#
# usage: gawk -v query_ind=<int> -v query_start=<int> -v ref_end=<int> -f find_pos_btop_minus.awk <<< <BTOP string>

{
    btop = $0
    i_query = query_start
    i_ref = ref_end
    
    while (length(btop) != 0) {
        regex_found = 0
        
        if (match(btop, /^([A-Z-])([A-Z-])/, arr)) {
            btop = substr(btop, RLENGTH + 1)
            regex_found = 1
            
            
            if (arr[2] == "-") {
                # insertion
                i_query += 1
                insertion_length += 1
                
                if (i_query == query_ind) {
                    # count adjacent insertions
                    match(btop, /^([^-]-)*/)
                    printf i_ref - insertion_length / (insertion_length + RLENGTH + 1)"\n"
                    exit
                }
            }
            else if (arr[1] == "-") {
                # deletion
                i_ref -= 1
                insertion_length = 0
            }
            else {
                # mismatch
                i_query += 1
                i_ref -= 1
                insertion_length = 0
                
                if (i_query == query_ind) {
                    printf i_ref"\n"
                    exit
                }
            }
        }
        
        if (match(btop, /^([0-9]+)/, arr)) {
            # match(es) (or padding)
            btop = substr(btop, RLENGTH + 1)
            regex_found = 1
            
            i_query += arr[1]
            i_ref -= arr[1]
            insertion_length = 0
            
            if (i_query >= query_ind) {
                printf i_ref - query_ind + i_query"\n"
                exit
            }
        }
        
        if (regex_found == 0) {
            print "Error: cannot parse BTOP string" > "/dev/stderr"
            printf "-1\n"
            exit 1
        }
    
    }
    print "Error: cannot calculate reference position" > "/dev/stderr"
    printf "-1\n"
    exit 1
}
