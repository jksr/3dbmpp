# 3DBMPP: 3D Beta-barrel Membrane Protein Predictor

### Overview:
Pipeline for 3D structure prediction for transmembrane beta barrel proteins (TMBs) or outer membrane proteins (OMPs).
The method was described in 
`High-resolution structure prediction of β-barrel membrane proteins.
W Tian, M Lin, K Tang, J Liang, H Naveed. 
*Proceedings of the National Academy of Sciences* 115 (7), 1511-1516`


### Note:
This software depends on rough information of secondary structure of a TMB.
This information can be obtained through 3rd-party software such as the ones listed in 
[http://www.ompdb.org/links.php](http://www.ompdb.org/links.php).

Due to the computation complexity, this package only provides the structure prediction for the barrel domains of the TMBs.
For the loop sampling mentioned in Tian *et. al*, please refer to
[https://github.com/uic-lianglab/ompg-public](https://github.com/uic-lianglab/ompg-public).

### Dependency:
This software depends the following python packages
* [numpy](https://numpy.org/)
* [biopython](https://biopython.org/)

and following software
* [BBQ algorithm](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20624) (BBQ algorithm is a part of [BioShell](http://www.bioshell.pl/bioshell.html))
* [Scwrl4](http://dunbrack.fccc.edu/scwrl4)

### Installation:

#### Linux:

`git clone https://github.com/jksr/3dbmpp-pipe.git`

`cd bin/src`

`make`

`cd ../..`

*Scwrl4* should also be installed. Please refer to [http://dunbrack.fccc.edu/scwrl4](http://dunbrack.fccc.edu/scwrl4). 


### Usage:

1. Fisrt, a folder needs to be created to store all the input files and results. 

2. Put the fasta file into the folder created.

3. Put a file with rough information of the starting and ending points of beta strands into the folder. Please see *example/1bxw.strands* for an example. The seqid of the starting and ending points shall be consistent with the fasta file. The staring and ending points of beta strands can be estimated via the secondary structure prediction tools listed in [http://www.ompdb.org/links.php](http://www.ompdb.org/links.php).

4. Information from the [*PSICOV*](https://github.com/psipred/psicov) sequence convariation analysis may help the structure prediction. Please refer to [https://github.com/psipred/psicov](https://github.com/psipred/psicov) for the installation and usage. This input is not obligatory. One can create a empty file with filename ending with .psicov in the folder to skip this step. However, the accuray of the prediction may be affected.

5. In our method, we classified TMBs into five groups.
One need to determine which group the pretoin is in before the predition. 
The five groups are listed below:

| Groups | Description                                      | Example PBD ids |
| ------ | ------------------------------------------------ | --------------- |
|  1     | Small TMBs (strand#<16) w/o inplugs or outclamps | 1bxw, 1qj8, 1p4t, 2f1t, 1thq, 2erv, 2lhf,  3dzm, 1qd6, 2f1c, 1k24, 1i78, 2wjr, 4pr7 |
|  2     | Small TMBs (strand#<16) w/ inplugs or outclamps  | 1t16, 1uyn, 1tly, 3aeh, 3bs0, 3dwo, 3fid, 3kvn, 4e1s |
|  3     | Medium oligomeric TMBs (16≤strand#<20)           | 2mpr, 1a0s, 2omf, 2por, 1prn, 1e54, 2o4v, 3vzt, 4n75 |
|  4     | Medium monomeric TMBs (16≤strand#<20)            | 2qdz, 2ynk, 3emn, 3rbh, 3syb, 3szv, 4c00, 4gey |
|  5     | Large TMBs (strand#≥20)                          | 1fep, 2fcp, 1kmo, 1nqe, 1xkw, 2vqi, 3csl, 3rfz, 3v8x, 4q35 |

6. Run the follow command to predict the 3D structure of the transmembrane beta barrel protein.

`python 3dbmpp.py --group *group_id* --folder *input_folder*  --scwrlpath "*path_to_scwrl4_executable*"`

For example, 

`python 3dbmpp.py --group 1 --folder example  --scwrlpath "*path_to_scwrl4_executable*"`

