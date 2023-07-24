write pseudo code for the python function below:
def _call_peaks_localmax(tx_id, count_signal, tx_reads, min_peak_length, max_peak_length, bin_width, min_cov, permutate_pval, poisson_pval, decay, mode, max_iter, seed, boundary,recursive_maxima_elimination):
    '''
    Input: bw signals,
    Return: local or global maxima position
    '''
    np.random.seed(seed)
    x = np.array(count_signal)
    gene_length = x.shape[0] 
    tx_global = [tx_id,0,gene_length,"+",tx_id]
    tx_global_reads = read_tid_frag_from_records(tx_global, tx_reads) 
    res = do_permutation(transcr=tx_global, read_transcript=tx_global_reads, max_iter=max_iter, pval_cutoff=permutate_pval, min_cov=min_cov, seed=seed)
    tx_cov, _height_to_pval = res[0:2]
    tx_cov = max(tx_cov, min_cov)
    half_bin_width = bin_width//2
    n_bins = max(1, gene_length//half_bin_width)
    bin_cov = np.zeros(n_bins)
    for i in range(n_bins):
        bin_cov[i] = np.mean(x[(i*half_bin_width):min((i + 2)*half_bin_width, gene_length)]) 
    cand_bins = np.nonzero(bin_cov > min_cov)[0] 
    cand_bins_list = [[i*half_bin_width,min((i + 2)*half_bin_width, gene_length)] for i in cand_bins] 
    cand_bins_merge = mergeIntervals(cand_bins_list) 
    all_maximas = []
    all_minimas = []
    for cand_bin_merge_index in range(len(cand_bins_merge)):
        this_bin_merge_start,this_bin_merge_end = cand_bins_merge[cand_bin_merge_index]
        this_maxima = find_local_maxima(arr=x[this_bin_merge_start:this_bin_merge_end]) 
        this_indexed_maxima = sorted([[i,int(val)] for i,val in enumerate(this_maxima) if val>0]) 
        this_indexed_maxima_filter = [i for i in this_indexed_maxima if i[1]>min_cov] 
        this_indexed_minima_filter = find_one_minima_between_maxima(arr=x[this_bin_merge_start:this_bin_merge_end], maxima_list=this_indexed_maxima_filter)
        last_num = len(this_indexed_maxima_filter)
        this_indexed_maxima_filter_merge = this_indexed_maxima_filter 
        this_indexed_minima_filter_merge = this_indexed_minima_filter 
        if last_num <=1:
            print("maxs_list size <=1, not iterable, not changed")
        elif not recursive_maxima_elimination:
            print("recursive_maxima_elimination == False, not changed")
        else:  
            last_num = 0 
            iter_num = 0
            while (len(this_indexed_maxima_filter_merge) != last_num) and (len(this_indexed_maxima_filter_merge) > 1): 
                last_num = len(this_indexed_maxima_filter_merge)
                this_indexed_maxima_filter_merge = merge_maximas_by_valley_cov(maxs_list=this_indexed_maxima_filter_merge,mins_list=this_indexed_minima_filter_merge,decay=decay) 
                this_indexed_minima_filter_merge = find_one_minima_from_multi_minima(maxima_list=this_indexed_maxima_filter_merge, mins_list=this_indexed_minima_filter) 
                iter_num += 1
        this_indexed_maxima_filter_merge = [[i[0]+this_bin_merge_start,i[1]]for i in this_indexed_maxima_filter_merge]
        try: 
            this_indexed_minima_filter_merge = [[i[0]+this_bin_merge_start,i[1]]for i in this_indexed_minima_filter_merge] 
        except: 
            this_indexed_minima_filter_merge = [[]]
        all_maximas.append(this_indexed_maxima_filter_merge) 
        all_minimas.append(this_indexed_minima_filter_merge)
    all_maximas_flat = [element for sublist in all_maximas for element in sublist]
    n_cand_peaks = len(all_maximas_flat)
    n_cand_peak = 0
    locals_list = [[0,0]] 
    fold_list = [] 
    peaks = []    
    for cand_bin_merge_index in range(len(cand_bins_merge)):
        this_bin_merge_start,this_bin_merge_end = cand_bins_merge[cand_bin_merge_index]
        tmp_maximas = all_maximas[cand_bin_merge_index]
        tmp_minimas = all_minimas[cand_bin_merge_index]
        left_bound = this_bin_merge_start
        left_valley = min_cov 
        for peak_index in range(len(tmp_maximas)): 
            try:
                right_bound,right_valley = tmp_minimas[peak_index] 
            except:
                right_bound = this_bin_merge_end
                right_valley = min_cov 
            try:
                max_index,local_max = tmp_maximas[peak_index]
            except:
                left_bound = right_bound
                left_valley = right_valley
                n_cand_peak += 1
                continue
            if max_index >= gene_length:
                max_index = gene_length-1
            if x[max_index] <= min_cov:
                left_bound = right_bound
                left_valley = right_valley
                n_cand_peak += 1
                continue
            if mode=="local":
                if max_index>locals_list[-1][0] and max_index<locals_list[-1][1]: 
                    s = locals_list[-1][0]
                    e = locals_list[-1][1]
                    fold = fold_list[-1]
                else:
                    most_left_start = max_index
                    while (most_left_start > left_bound) and (x[most_left_start - 1] >= min_cov): 
                        most_left_start -= 1
                    most_right_end = max_index
                    while (most_right_end < (gene_length - 1)) and (x[most_right_end + 1] >= min_cov): 
                        most_right_end += 1
                    locals_list.append([most_left_start,most_right_end]) 
                    subset = min(500,(most_right_end-most_left_start)) 
                    gini_val = np.nan_to_num(gini(np.random.choice(x[most_left_start:most_right_end], subset, replace=False)))
                    print(datetime.datetime.now(),"\t\t","gini: ",gini_val) 
                    fold = int(3+47*(1-math.log2(gini_val+1))) 
                    fold_list.append(fold)
                    span = fold*(most_right_end-most_left_start)
                    if most_right_end-most_left_start <= 10:
                        print(datetime.datetime.now(),"\t\t","initial local <= 10, skip this bin !")
                        left_bound = right_bound
                        left_valley = right_valley
                        n_cand_peak += 1
                        continue
                    s = int(0.5*(most_left_start+most_right_end)-0.5*span)
                    e = int(s+span)
                    if gene_length <= span:
                        s = 0
                        e = gene_length
                    else:
                        if s < 0: 
                            s = 0 
                            e = span
                        elif e > gene_length:
                            s = gene_length - span
                            e = gene_length
                    tx_local = [tx_id,s,e,"+",tx_id]
                    tx_local2 = [tx_id,most_left_start,most_right_end,"+",tx_id]
                    tx_local_reads = read_tid_frag_from_records(tx_local2, tx_reads) 
                    res = do_permutation(transcr=tx_local, read_transcript=tx_local_reads, max_iter=max_iter, pval_cutoff=permutate_pval, min_cov=min_cov, seed=seed) 
                    tx_cov, _height_to_pval = res[0:2]
                    tx_cov = max(tx_cov, min_cov)
            elif mode=="global":
                s = 0
                e = gene_length 
            else:
                raise Exception("only global and local mode for background cov estimation supported !")
            if local_max-tx_cov <= 0:
                left_bound = right_bound
                left_valley = right_valley
                n_cand_peak += 1
                continue
            start = end = max_index   
            max_index_new = max_index
            if boundary == "localmaxdecay":
                while (start > left_bound) and (x[start - 1] >= decay*local_max) and (x[start - 1] >= left_valley): 
                    start -= 1
                while (end < right_bound-1) and (x[end + 1] >= decay*local_max) and (x[end + 1] >= right_valley): 
                    end += 1
            elif boundary == "background": 
                while (start > left_bound) and (x[start - 1] >= tx_cov) and (x[start - 1] >= left_valley): 
                    start -= 1
                while (end < right_bound-1) and (x[end + 1] >= tx_cov) and (x[end + 1] >= right_valley): 
                    end += 1
            else:
                raise Exception("only localmaxdecay and background mode for boundary are supported !")
            tx_peak = [tx_id,start,end,"+",tx_id]
            tx_peak_reads = read_tid_frag_from_records(tx_peak, tx_reads)
            if mode=="local":
                fold = max(10,fold) 
                half_peak_len = int((fold-1)*0.5*(end-start)) 
                tmp_local_start = max(0,start-half_peak_len)
                tmp_local_end = min(end+half_peak_len,gene_length)
                tx_peak_local = [tx_id,tmp_local_start,tmp_local_end,"+",tx_id]
                tx_peak_local_reads = read_tid_frag_from_records(tx_peak_local, tx_reads)
                poisson_p = poissonP(reads_in_gene=len(tx_peak_local_reads), reads_in_peak=len(tx_peak_reads), gene_length=tmp_local_end-tmp_local_start, peak_length=end-start)
            elif mode=="global":
                poisson_p = poissonP(reads_in_gene=len(tx_global_reads), reads_in_peak=len(tx_peak_reads), gene_length=gene_length, peak_length=end-start)
            else:
                raise Exception("only global and local mode supported !")
            print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"poisson_p:",poisson_p)
            if ((end - start) < min_peak_length) or ((end - start) > max_peak_length): 
                print(datetime.datetime.now(),"\t\t","finished (lenth",(end - start),"not pass filter)",n_cand_peak,"/",n_cand_peaks,"cand_bins !")
            elif (poisson_p > poisson_pval):
                print(datetime.datetime.now(),"\t\t","finished (poisson_pval",poisson_p,"not pass filter",poisson_pval,")",n_cand_peak,"/",n_cand_peaks,"cand_bins !")
            else:
                if (len(_height_to_pval) >= 1):
                        try:
                            tmp_fdr = _height_to_pval[int(local_max)]
                        except:
                            if int(local_max) >= max(_height_to_pval.keys()):
                                tmp_fdr = float(0.0)
                            else:
                                tmp_fdr = float('nan')
                else:
                    tmp_fdr = float('nan')
                peaks.append([tx_id, start, end, local_max, max_index_new, tx_cov, tmp_fdr, poisson_p]) 
            left_bound = right_bound
            left_valley = right_valley
            n_cand_peak += 1
    return peaks

