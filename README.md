# AmproCode-scripts

A database search algorithm proving the theoretical feasiblility of AmproCode and computational simulation to estimate the coverage of the whole proteome and secretome

## Usage

### Searching from fasta

All scripts are located in "scripts", the first option is the sequence database(fasta), followed by "code" of the sample. For example:

`python search_CKMDE.py ../databases/secreted_seq.fasta 0.97 1 0.98 0.49`

> rank\= 1 P0DMC3|ELA\_HUMAN 7.468776545216382e-05 1
> 
> rank\= 2 Q9BYW3|DB126\_HUMAN 0.02542741131229176 1
> 
> rank\= 3 P16860|ANFB\_HUMAN 0.029564160560074715 1
> 
> rank\= 4 P10092|CALCB\_HUMAN 0.038327614363006135 1
> 
> rank\= 5 P11686|PSPC\_HUMAN 0.05956680608347864 1
>
> ... ...

`python search_CKMDEY.py ../databases/secreted_seq.fasta 2.06 0.01 1.0 0.95 1.00`

> rank\= 1 P10997|IAPP\_HUMAN 0.000384893975971079 2
> 
> rank\= 2 Q765I0|UTS2B\_HUMAN 0.000384893975971079 2
> 
> rank\= 3 Q6ZRU5|YQ032\_HUMAN 0.011808715519036728 1
> 
> rank\= 4 P60022|DEFB1\_HUMAN 0.01356380669196311 1
> 
> rank\= 5 P04808|REL1\_HUMAN 0.014103506613592942 1
>
> ... ...

`python search_CKMDE.py ../databases/UP000005640_9606.fasta 6 2 2 11`

> rank\= 1 sp|P59666 0 2
> 
> rank\= 2 sp|P01308 0 2
> 
> rank\= 3 sp|Q4KMG9 0.00046825665561844865 1
> 
> rank\= 4 sp|Q9BY78 0.0006505979351558722 1
> 
> rank\= 5 sp|Q76LX8 0.000841947322775427 1
>
> ... ...

Each line of the ouput: rank, protein\_name, **cos\_distance**, degeneracy

Note: cos\_similarity \= 1 - cos\_distance

### Simulation for noises

Adding random noise for “code” of each protein in proteome, the first option is the sequence database(fasta), followed by sigma of the noise(draw from gaussian). For example:

`python comp2seq_CK.py ../databases/secreted_seq.fasta 0.01`&#x20;

> 2578 top1\= 0.06361520558572537 top3\= 0.15438324282389448

Output: number of code, probabilityof the correct answer appearing in the top1 ranking, probabilityof the correct answer appearing in the top3 ranking.
