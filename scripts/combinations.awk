#https://stackoverflow.com/questions/56911991/how-would-i-loop-over-pairs-of-values-without-repetition-in-bash
###################
# Calculate all combinations of a set of strings, see
# https://rosettacode.org/wiki/Combinations#AWK
###################
# Usage:
#$ awk -f combinations.awk <<< '1 2 3 4'
#1 2
#1 3
#1 4
#2 3
#2 4
#3 4

function get_combs(A,B, i,n,comb) {
    ## Default value for r is to choose 2 from pool of all elements in A.
    ## Can alternatively be set on the command line:-
    ##    awk -v r=<number of items being chosen> -f <scriptname>
    n = length(A)
    if (r=="") r = 2

    comb = ""
    for (i=1; i <= r; i++) { ## First combination of items:
        indices[i] = i
        comb = (i>1 ? comb OFS : "") A[indices[i]]
    }
    B[comb]

    ## While 1st item is less than its maximum permitted value...
    while (indices[1] < n - r + 1) {
        ## loop backwards through all items in the previous
        ## combination of items until an item is found that is
        ## less than its maximum permitted value:
        for (i = r; i >= 1; i--) {
            ## If the equivalently positioned item in the
            ## previous combination of items is less than its
            ## maximum permitted value...
            if (indices[i] < n - r + i) {
                ## increment the current item by 1:
                indices[i]++
                ## Save the current position-index for use
                ## outside this "for" loop:
                p = i
                break
            }
        }

        ## Put consecutive numbers in the remainder of the array,
        ## counting up from position-index p.
        for (i = p + 1; i <= r; i++) indices[i] = indices[i - 1] + 1

        ## Print the current combination of items:
        comb = ""
        for (i=1; i <= r; i++) {
            comb = (i>1 ? comb OFS : "") A[indices[i]]
        }
        B[comb]
    }
}

# Input should be a list of strings
{
    split($0,A)
    delete B
    get_combs(A,B)
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (comb in B) {
        print comb
    }
}
