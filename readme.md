- [eXtended Contact Mode Score (xCMS)](#sec-1)
  - [Methods](#sec-1-1)
  - [Pipeline](#sec-1-2)
  - [xCms back compatible with Cms](#sec-1-3)
  - [Statistical significance of cms](#sec-1-4)
  - [Statistical significance of xcms value](#sec-1-5)
    - [Run on a random dataset](#sec-1-5-1)
- [Assess AutoVina on DUDE](#sec-2)
  - [Presumption](#sec-2-1)
    - [Prove](#sec-2-1-1)
  - [Use mutiple template pockets and ligands](#sec-2-2)
  - [Random prediction](#sec-2-3)

# eXtended Contact Mode Score (xCMS)<a id="orgheadline7"></a>

Extend the application of Contact Mode Score to non-identical complexes

## Methods<a id="orgheadline1"></a>

1.  Contact Matrix
    1.  Corresponding residues of two binding pockets
        1.  [X] <a id="orgradiotarget1">subject</a>
            1.  [X] Ligands of two pockets have a Tc >= 0.5
            2.  [X] Sequence similarity < 30% between the two associated proteins
            3.  [X] shaore >= 50 atmic ligand-protein contacts of same type
            4.  [X] the pdb files provided by Apoc contains only C\_alpha atoms
    2.  [X] coordinates of the effective points
        1.  use the pdb files downloaded from PDB to grep the proteins' residues
        2.  use the mixed resolution as in the CMS
        3.  use the <a id="orgradiotarget2">virtual effective</a> points to compensate for miss-matched residues
            1.  the contact between a ligand atom and a virtual effective point is always **None**
        4.  [X] retrieve the ligand from mmCIF file and convert to mol2 format
            1.  [X] can tag such as "3zyc\_GCP\_A\_1749" uniquely identify a ligand from the pdb file?
                Yes, but it is not safe to do so in a pdb file
            2.  [X] check out
                1.  [X] Python PDBx
                2.  [X] BioPython
                3.  [X] OpenBabel c++ library
                    1.  residue names become "LIG" when read by babel
                    2.  [X] how babel parse the mmCIF file?
                        hacked the parser in babel for mmCIF in ordr to retrieve auth\* fields
    3.  [ ] certain ligands cannot be grepped
        e.g. 3q2k\_NAD\_P\_500
    4.  [ ] certain pdbs miss some atoms
        1.  in /ddnB/work/jaydy/dat/pdb/vn/pdb3vn9.ent
            
            <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
            
            
            <colgroup>
            <col  class="org-left" />
            
            <col  class="org-right" />
            
            <col  class="org-left" />
            
            <col  class="org-left" />
            
            <col  class="org-left" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-left" />
            </colgroup>
            <tbody>
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">674</td>
            <td class="org-left">N</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">130</td>
            <td class="org-right">-2.316</td>
            <td class="org-right">-44.226</td>
            <td class="org-right">2.992</td>
            <td class="org-right">1.00</td>
            <td class="org-right">90.63</td>
            <td class="org-left">N</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">675</td>
            <td class="org-left">CA</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">130</td>
            <td class="org-right">-1.581</td>
            <td class="org-right">-44.192</td>
            <td class="org-right">4.257</td>
            <td class="org-right">1.00</td>
            <td class="org-right">83.49</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">676</td>
            <td class="org-left">C</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">130</td>
            <td class="org-right">-1.888</td>
            <td class="org-right">-42.913</td>
            <td class="org-right">5.017</td>
            <td class="org-right">1.00</td>
            <td class="org-right">100.53</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">677</td>
            <td class="org-left">O</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">130</td>
            <td class="org-right">-2.031</td>
            <td class="org-right">-41.844</td>
            <td class="org-right">4.392</td>
            <td class="org-right">1.00</td>
            <td class="org-right">97.24</td>
            <td class="org-left">O</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">678</td>
            <td class="org-left">CB</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">130</td>
            <td class="org-right">-0.072</td>
            <td class="org-right">-44.310</td>
            <td class="org-right">4.027</td>
            <td class="org-right">1.00</td>
            <td class="org-right">73.44</td>
            <td class="org-left">C</td>
            </tr>
            </tbody>
            </table>
        2.  in /ddnB/work/jaydy/dat/pdb/t9/pdb3t9a.ent
            
            <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
            
            
            <colgroup>
            <col  class="org-left" />
            
            <col  class="org-right" />
            
            <col  class="org-left" />
            
            <col  class="org-left" />
            
            <col  class="org-left" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-right" />
            
            <col  class="org-left" />
            </colgroup>
            <tbody>
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1612</td>
            <td class="org-left">N</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">8.315</td>
            <td class="org-right">-2.250</td>
            <td class="org-right">6.872</td>
            <td class="org-right">1.00</td>
            <td class="org-right">21.85</td>
            <td class="org-left">N</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1613</td>
            <td class="org-left">CA</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">9.731</td>
            <td class="org-right">-2.356</td>
            <td class="org-right">7.160</td>
            <td class="org-right">1.00</td>
            <td class="org-right">23.40</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1614</td>
            <td class="org-left">C</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">10.257</td>
            <td class="org-right">-0.984</td>
            <td class="org-right">7.551</td>
            <td class="org-right">1.00</td>
            <td class="org-right">23.08</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1615</td>
            <td class="org-left">O</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">9.910</td>
            <td class="org-right">0.026</td>
            <td class="org-right">6.922</td>
            <td class="org-right">1.00</td>
            <td class="org-right">22.77</td>
            <td class="org-left">O</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1616</td>
            <td class="org-left">CB</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">10.506</td>
            <td class="org-right">-2.904</td>
            <td class="org-right">5.947</td>
            <td class="org-right">1.00</td>
            <td class="org-right">23.90</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1617</td>
            <td class="org-left">CG</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">10.539</td>
            <td class="org-right">-1.968</td>
            <td class="org-right">4.735</td>
            <td class="org-right">1.00</td>
            <td class="org-right">27.61</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1618</td>
            <td class="org-left">CD</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">11.471</td>
            <td class="org-right">-2.450</td>
            <td class="org-right">3.608</td>
            <td class="org-right">1.00</td>
            <td class="org-right">32.44</td>
            <td class="org-left">C</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1619</td>
            <td class="org-left">OE1</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">12.400</td>
            <td class="org-right">-3.248</td>
            <td class="org-right">3.888</td>
            <td class="org-right">1.00</td>
            <td class="org-right">35.07</td>
            <td class="org-left">O</td>
            </tr>
            
            
            <tr>
            <td class="org-left">ATOM</td>
            <td class="org-right">1620</td>
            <td class="org-left">OE2</td>
            <td class="org-left">GLU</td>
            <td class="org-left">A</td>
            <td class="org-right">238</td>
            <td class="org-right">11.268</td>
            <td class="org-right">-2.020</td>
            <td class="org-right">2.448</td>
            <td class="org-right">1.00</td>
            <td class="org-right">32.92</td>
            <td class="org-left">O</td>
            </tr>
            </tbody>
            </table>
        3.  problem resolved if by using virtual effective points
    5.  list file from apoc contains weird keys words such as "### XYM 35"
        1.  original
            /ddnB/work/jaydy/dat/apoc/RS2.lst
        2.  cleaned
            /work/jaydy/dat/apoc/RS2\_valid.lst
    6.  Corresponding heavy atoms of two ligands
        1.  Kcombu (currently used)
        2.  SIMCOMP (used in Apoc paper)
            <a id="orgradiotarget3">SIMCOMP</a> tend to yield higher Tc than <a id="orgradiotarget4">pkcombu</a>

## Pipeline<a id="orgheadline2"></a>

1.  Original pdb file for the complex
    1.  apoc\_inputs.PdbPathTask
    2.  apoc\_inputs.DecompressedPdb
2.  Grep Ligand atoms from the pdb
    1.  apoc\_inputs.LigandExpStructureInPdb
3.  Run apoc
    1.  run\_control\_apoc.LpcApocResultTask
4.  Run kcombu
    1.  run\_control\_apoc.LpcKcombuResult
5.  Parsers
    1.  run\_control\_apoc.ApocResultParer
    2.  run\_control\_apoc.PkcombuAtomMatchParser
6.  Run xcms
    1.  xcms.LpcApocXcms
7.  Collect data for analysis
    1.  collect\_xcms.AtomicXcmsCollection
    2.  collect\_xcms.AtomicXcmsTable

## xCms back compatible with Cms<a id="orgheadline3"></a>

The value of xCms should be <a id="orgradiotarget5">back-compatible</a> with Cms.

1.  materials
    complexes in the Astext dataset
2.  methods
    1.  apply rotation and translation to each ligand
        1.  to calculate Cms for each complex
            1.  calculate the Cms between the variational conformations and the native conformation
        2.  to calculate xCms for eac complex
            1.  find the neighboring residues in contact with ligand atoms
            2.  align the neighboring residues using Apoc
            3.  for the matching residues from Apoc, add virtual points
            4.  calculate the xCms between the variational conformations and the native conformation

## Statistical significance of cms<a id="orgheadline4"></a>

1.  run MC for each complex in astex with weak constraints
    1.  vdw \* 0.5
2.  cluster the conformations based on their translational and rotational vector
    1.  [X] experiment with DBSCAN
        1.  do not use DBSCAN becuase it removes the outliers
    2.  [X] experiment with KMeans
        1.  Kmeans maintain the density distribution as the original ponts
    3.  experiment with the grid method
        1.  [X] grid method provides evenly distributed points
        2.  [ ] in bining the points, use the median point as the center
            1.  Be careful with rotational angles; they **DO NOT** conform to vector calculation
                1.  [ ] experiment with not bining the rotational angles
3.  calculate pair-wise p-value
4.  sample
    1.  [X] stratefied sampling
        Assume that the # of grid cells is n,
        the # of samples from each complex is &prop; to n^2
5.  [ ] dependence on the ligand size and the protein size
6.  [ ] analyze the distribution of the cms values

## Statistical significance of xcms value<a id="orgheadline6"></a>

### Run on a random dataset<a id="orgheadline5"></a>

1.  Random set (<a id="orgradiotarget6">RS2</a>) from Apoc
    /ddnB/work/jaydy/dat/apoc/RS2.lst
    /ddnB/work/jaydy/dat/apoc/RS2\_valid.lst
2.  script
    ./src/run\_lpc\_randomset.py
3.  [X] run time errors
    Cannot do much about the errors below because of the in-perfection of the source of data.
    1.  missing pdb files
        /ddnB/work/jaydy/dat/apoc/missing\_pdbs.txt
    2.  protein pairs that have mis-matched elements
        /ddnB/work/jaydy/dat/apoc/mis\_matched\_prt.txt
    3.  kcombu failures
        /ddnB/work/jaydy/dat/apoc/kcombu\_failures.txt
4.  [ ] curate <a id="orgradiotarget7">dataset</a>
    1.  [ ] calculate the all-against-all apoc, kcombu using the successful pockets dataset
    2.  [ ] cluster using DBSCAN based on the ps-score and Tc between two pockets
        1.  \sqrt{ps-score^2 + Tc^2} &#x2013;> similarity
        2.  or use the product directly
        3.  similarity(ps-score = 0.4, Tc = 0.4) &#x2013;> eps
    3.  [ ] use the centroids' xcms to calculate the statistical significance
5.  curate subject
    ./src/curate.py
    ./src/run\_curate.py
    1.  dbscan for subject dataset
        
        <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
        
        
        <colgroup>
        <col  class="org-right" />
        
        <col  class="org-right" />
        </colgroup>
        <thead>
        <tr>
        <th scope="col" class="org-right">eps</th>
        <th scope="col" class="org-right">#clusters</th>
        </tr>
        </thead>
        
        <tbody>
        <tr>
        <td class="org-right">1.0</td>
        <td class="org-right">1070</td>
        </tr>
        
        
        <tr>
        <td class="org-right">1.1</td>
        <td class="org-right">764</td>
        </tr>
        
        
        <tr>
        <td class="org-right">1.15</td>
        <td class="org-right">620</td>
        </tr>
        
        
        <tr>
        <td class="org-right">1.2</td>
        <td class="org-right">406</td>
        </tr>
        
        
        <tr>
        <td class="org-right">1.3</td>
        <td class="org-right">103</td>
        </tr>
        </tbody>
        </table>

# TODO Assess AutoVina on DUDE<a id="orgheadline12"></a>

## DONE Presumption<a id="orgheadline9"></a>

Similar ligand and pockets yield high spearmanr between the distance vectors

### DONE Prove<a id="orgheadline8"></a>

1.  [X] Use BioLip's structure to query
2.  [X] Cluster the ligands in BioLip
3.  [X] do the clustering of the proteins using CD-HIT?
    <http://weizhongli-lab.org/cd-hit/>
    
    ```python
    # generate fasta sequence from pdb file
    import pybel
    prt = pybel.readfile('pdb', "/ddnB/work/jaydy/dat/BioLip/prt/01/101mA.pdb").next()
    fasta_str = prt.write('fasta')
    ```
    
    No, because we will rank the templates
4.  [X] find how spearmanr evolves with ps\_score and TC

## DONE Use mutiple template pockets and ligands<a id="orgheadline10"></a>

1.  Non-redundant dataset of proteins and ligands
    <http://zhanglab.ccmb.med.umich.edu/BioLiP/download.html>
    1.  filter the ligands with
        1.  size constraints (6 < size)
        2.  [ ] tanimoto coefficients > 0.5 (using FastSearch in OpenBabel <http://openbabel.org/wiki/FastSearch>)
            Note that pkcombu and babel gives different tanimoto coefficients for the same ligand pair
    2.  [X] Index the ligands
        1.  [X] use only small ligands with #atom < 999
        2.  [X] /home/jaydy/Workspace/Bitbucket/xcms/src/ligand\_nr.bash
        3.  [X] /home/jaydy/Workspace/Bitbucket/xcms/src/apc\_pbs/wq\_ligand\_nr.pbs
    3.  [X] Index
        
        ```sh
        cd /ddnB/work/jaydy/dat/BioLip/
        babel ligand_nr.sdf -ofs
        # results are written to /ddnB/work/jaydy/dat/BioLip/ligand_nr.fs
        ```
    4.  [X] Search
        
        ```sh
        # example command
        cd /work/jaydy/dat/BioLip
        babel ligand_nr.fs outfile.sdf -s ./sml_ligand_nr/08/208l_CYS_A_1.pdb -at0.5
        ```
    5.  BioLip-referenced SpearmanR

## TODO Random prediction<a id="orgheadline11"></a>

AutoVina has an option for random prediciton of ligand within the binding site.
