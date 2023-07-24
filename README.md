<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
<!-- [![GPL-2.0 License][license-shield]][license-url] -->



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://gitee.com/hunnngry/cfPeak">
    <img src="_plots/cfPeak.png" alt="Logo" width="160" height="80">
  </a>

<h2 align="center">cfPeak</h2>

  <p align="center">
    cfPeak, peak analysis finds recurrently protected narrow regions with clinical potential
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#References">References</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->

## About The Project
we report cfPeak, a new computational method and data pipeline with multiple advantages for peak analysis in cell-free RNA-seq, which is also potentially applicable to other data types, including traditional cellular CLIP-seq or RIP-seq.

<!-- PROJECT pipeline -->
<br />
<p align="center">
  <a href="https://gitee.com/hunnngry/cfPeak">
    <img src="_plots/pipeline.png" alt="Pipeline" width="480" height="400">
  </a>


<!-- GETTING STARTED -->

## Getting Started

Following steps below to get a local copy and run

### Prerequisites
* slurm or lsf cluster system
* PC or server with >=8 cores and >=16 GB memory
* anaconda or miniconda installed (https://www.anaconda.com/download)

### Installation

1. Clone the repo
   ```sh
   # from gitee
   git clone https://gitee.com/hunnngry/cfPeak.git
   # or from github
   git clone https://github.com/hunnngry/cfPeak.git
   ```
2. Create conda global environment
   ```sh
   cd cfPeak
   # The default conda solver is a bit slow and sometimes has issues with selecting the special version packages. We recommend to install mamba as a drop-in replacement
   conda install -c conda-forge mamba
   mamba env create -n cfpeak -f ./snakemake/envs/cfpeak.yml
   ```
3. Create snakemake pipeline environment

   a. required for all
   ```sh
   # activate global environment
   source activate cfpeak
   # if you need to run cell-free smRNA-seq:
   snakemake --use-conda --conda-create-envs-only -j 1 --configfile config/test_small.yaml -s snakemake/call_peaks.snakemake
   # if you need to run cell-free totalRNA-seq:
   snakemake --use-conda --conda-create-envs-only -j 1 --configfile config/test_long.yaml -s snakemake/call_peaks_long.snakemake
   # some other config and snakemake files are also provided to meet different input types 
   ```

   b. optional: if you need to run Piranha (v1.2.1)
   * By default, basical env (snakemake/envs/piranha.yml) that needed to run Piranha has already installed in last step, but need install countreg package manually to make it run: this R package unfortunately is not yet on conda, so must be installed manually in a sort of hacky way. Locate your installed conda environments (by default in .snakemake/conda/[hash], as described in Snakemake documentation). Find out which [hash].yml file corresponds with the environment named piranha. Load that environment with conda activate .snakemake/conda/[hash]. Then run R to enter the R shell and install.packages("countreg", repos="http://R-Forge.R-project.org") to install countreg.
   * We also provide two optional piranha version: in-house adapted cell-free Piranha software and simplified Piranha pkg in R (default). Users could switch to the former Piranha version in peak_common.snakemake by commenting rule call_peaks_piranha and uncommenting rule call_peaks_piranha2 if needed.

   c. optional: if you need to run CLIPper (v2.1.2) 
   * install CLIPper according to https://github.com/YeoLab/clipper
   * modify to suit tx mode:
     * find installed path of CLIPper call peak python script (e.g., ~/anaconda3/envs/clipper3/lib/python3.7/site-packages/clipper/src/call_peak.py)
     * mask 927-929 rows (or else it will add "chr" prefix in seqnames and cause error)
     * find main python script (e.g., ~/anaconda3/envs/clipper3/lib/python3.7/site-packages/clipper/src/main.py)
     * remove  comment of  203,204,209 rows (parameters: reverse_strand False  max_width  min_width )
     * add comment 232-234 rows
   * use hg38_tx reference provided in cfpeak (hg38txNoDNAnewTxID.gtf in https://cloud.tsinghua.edu.cn/f/9d8cee33da6e4aacbc40/?dl=1)
   * or you can add your own reference according to https://github.com/YeoLab/clipper/wiki/Supporting-additional-species
  
   d. optional: if you need to run CLAM (v1.2.0) 
   * install CLIPper according to https://github.com/YeoLab/clipper
   * modify to suit tx mode:
     * modify permutation_callpeak.py： row#439，(this fixed error for not considering left boundary less than 0 when defining flag: --extend)
     * permutation_peakcaller: count_pileup_heights(tlen, reads) row 307: # seem forgot to decision, here seem not center as tag position (this lead to a shift from real peak region)

4. Prepare/Download necessary genome/transcriptome reference file and annotation 
   * human hg38 chrom_sizes and annotation track files used in the article can be downloaded from https://cloud.tsinghua.edu.cn/f/9d8cee33da6e4aacbc40/?dl=1
   * or you can create your own reference of annotation file 


<!-- USAGE EXAMPLES -->

## Usage

1. Activate the created global environment
   ```sh
   source activate cfpeak
   export PATH=~/anaconda3/envs/cfpeak/bin:$PATH 
   dst="test_small"
   ```
2. Modify the config file based on your purpose
   detailed description are shown in config/test_small.yaml file
3.1 Run in server
   ```sh
    snakemake \
      --rerun-incomplete --keep-going --printshellcmds --reason --use-conda --nolock --latency-wait 20 --restart-times 1 --jobs 14 \
      --snakefile snakemake/call_peaks.snakemake \
      --configfile config/${dst}.yaml \
      > logs/${dst}/run-${dst}.log 2>&1 &
   ```
3.2 Run in cluster
   ```sh
    snakemake \
      --rerun-incomplete --keep-going --printshellcmds --reason --use-conda --nolock --latency-wait 20 --restart-times 1 --jobs 14 \
      --snakefile snakemake/call_peaks.snakemake \
      --configfile config/${dst}.yaml \
      --cluster-config snakemake/cluster_slurm.json \
      --cluster "sbatch --cpus-per-task 1 -n {cluster.threads} -J {cluster.jobname} -p {cluster.partition}  {cluster.resources} -o {cluster.output} -e {cluster.error}" \
      > logs/${dst}/run-${dst}.log 2>&1 &
   ```
4. Check output
   * cfPeak sample peak file: output/test_small/call_peak_all/cfpeak_by_sample/b5_d50_p1/s1.bed
   * cfPeak sample peak file (filtered by CNN): output/test_small/call_peak_all/cfpeakCNN_by_sample/b5_d50_p1/s1.bed
   * cfPeak consensus peak file: output/test_small/call_peak_all/cfpeak/b5_d50_p1.bed
   * cfPeak consensus peak count matrix: output/test_small/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1.txt


<!-- LICENSE -->

## License

Distributed under the GPL-2.0 License License. See `LICENSE` for more information.



<!-- CONTACT -->

## Contact

@hunnngry (bpf19@mails.tsinghua.edu.cn)

lulab (lulab1@tsinghua.edu.cn)


<!-- References -->

## References
* submitted, 2023, CfRNA peak analysis finds recurrently protected narrow regions with clinical potential.

<!-- ACKNOWLEDGEMENTS -->

## Acknowledgements
This work is supported by Tsinghua University Spring Breeze Fund (2021Z99CFY022), National Natural Science Foundation of China (81972798, 32170671, 81902384), National Key Research and Development Plan of China (2019YFC1315700), National Science and Technology Major Project of China (2018ZX10723204, 2018ZX10302205), Tsinghua University Guoqiang Institute Grant (2021GQG1020), Tsinghua University Initiative Scientific Research Program of Precision Medicine (2022ZLA003), Bioinformatics Platform of National Center for Protein Sciences (Beijing) (2021-NCPSB-005). This study was also supported by Beijing Advanced Innovation Center for Structural Biology, Bio-Computing Platform of Tsinghua University Branch of China National Center for Protein Sciences, Interdisciplinary Clinical Research Project of Peking University First Hospital and the Capital Health Research and Development of Special, Open Research Fund Program of Beijing National Research Center for Information Science and Technology.


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[contributors-shield]: https://img.shields.io/github/contributors/hunnngry/cfPeak.svg?style=for-the-badge

[contributors-url]: https://github.com/hunnngry/cfPeak/graphs/contributors

[forks-shield]: https://img.shields.io/github/forks/hunnngry/cfPeak.svg?style=for-the-badge

[forks-url]: https://github.com/hunnngry/cfPeak/network/members

[stars-shield]: https://img.shields.io/github/stars/hunnngry/cfPeak.svg?style=for-the-badge

[stars-url]: https://github.com/hunnngry/cfPeak/stargazers

[issues-shield]: https://img.shields.io/github/issues/hunnngry/cfPeak.svg?style=for-the-badge

[issues-url]: https://github.com/hunnngry/cfPeak/issues

<!-- [license-shield]: https://img.shields.io/github/license/hunnngry/cfPeak.svg?style=for-the-badge

[license-url]: https://github.com/hunnngry/cfPeak/blob/master/LICENSE.txt -->