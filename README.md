# PanGraphViewer -- show pangenome graphs in an easy way

![PyPI - Python Version](https://img.shields.io/pypi/pyversions/Django?color=green&logoColor=green&style=social) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/tf-chan-lab/panGraphViewer?color=green&logoColor=green&style=social) ![GitHub](https://img.shields.io/github/license/tf-chan-lab/panGraphViewer?color=green&logoColor=green&style=social) ![GitHub all releases](https://img.shields.io/github/downloads/tf-chan-lab/panGraphViewer/total?color=green&logoColor=green&style=social) 

### Table of Contents
<img align="right" width="300" src="doc/src/demo.jpg"/>

- [Versions and dependencies](#versions-and-dependencies)
  - [Desktop-based panGraphViewer](#desktop-based-pangraphviewer)
    - [Library installation for the desktop-based version](#library-installation-for-the-desktop-based-version)
    - [Start the desktop-based version](#start-the-desktop-based-version)
  - [Web-based panGraphViewer](#web-based-pangraphviewer)
    - [Library installation for the web-based version](#library-installation-for-the-web-based-version)
    - [Start the web-based version](#start-the-web-based-version)
- [Files accepted in the application](#files-accepted-in-the-application)
  - [GFA files](#gfa-files)
  - [VCF files](#vcf-file)
  - [Annotation files](#annotation-files)
- [Q&A:](#qa)
  - [Operating system](#operating-system)
  - [Minimum computational resource needed](#minimum-computational-resource-needed)
  - [Which application should I use](#which-application-should-i-use)
  - [How to use the application](#how-to-use-the-application)
  - [Backbone sample](#backbone-sample)
  - [Colors shown in the graph](#colors-shown-in-the-graph)
  - [Type of graphs](#type-of-graphs)
  - [Shapes shown in the graph](#shapes-shown-in-the-graph)
  - [Different variations](#different-variations)
  - [Pros and cons of the application](#pros-and-cons-of-the-application)

**Please read the information below carefully to have a better experience in using PanGraphViewer.**

---
## Versions and dependencies
Here we provide **two** application versions:
```
● Desktop-based application
● Web browser-based application
```
Overall, **Python3** is needed to run this software and we highly recommend using ``miniconda3`` to install all ``python3`` libraries. You may download ``miniconda3`` from the link as shown below
```
● On Windows systems, you can download miniconda3 at 
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

● On macOS systems, you can download miniconda3 at 
  https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh  

● On Linux systems,  you can download miniconda3 at 
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

After the installation of ``miniconda3``, you can follow the steps below to ensure `panGraphViewer` can be executed.

---
### Desktop-based panGraphViewer

#### Library installation for the desktop-based version

***Steps on different systems***

* If you use ``Windows`` systems, you may need to find or search ``Anaconda Powershell Prompt (miniconda3)`` first and then open it.

* If you use ``macOS`` or ``Linux`` system, you may open ``Terminal`` first and then type the command line below if you haven't set ``conda`` into your ``PATH``
    ```
    $ export PATH=/full/path/to/miniconda3/bin:$PATH  ## Please modify the path based on your ENV
    ```

After the steps above, you can install the ``python3`` libraries by typing:

```
conda config --add channels conda-forge
conda config --add channels bioconda
conda install pyqt pyqtwebengine configparser pandas bokeh==2.2.3 dna_features_viewer natsort attrdict networkx 
```

**Note:** 

1. On ``Linux`` or ``macOS`` system, ``pysam`` is needed. You may install this package using

    ```
    $ conda install pysam 
    ```

2. On ``Windows`` platforms, as `pysam` is not available, we use a windows-version ``samtools`` package instead. Additional libraries below are needed and can be installed using

    ```
    > conda install m2-base pyfaidx
    ```

3. As ``pip`` may introduce library conflicts for the ``desktop-based application``, we **do not recommend** using ``pip`` to install the needed python3 libraries, particularly for ``PyQtWebEngine``. However, we still provide such an option if you prefer to do so. 
   
   Firstly, you need to go to the ``panGraphViewerApp`` directory via ``Anaconda Prompt (miniconda3) Powershell`` or ``Terminal`` and then you may follow the command line below to install the python libraries
   ```
    pip install -r requirements.txt           ## On Linux or macOS systems
    pip install -r requirements_windows.txt   ## On Windows systems
   ```
    When you confront ``PyQtWebEngine`` calling problem after using ``pip`` to install, you may uninstall it and try to use ``conda`` to install. If the problem is still there, we are afraid that ``PyQtWebEngine`` is not compatible with your OS system. However, you can still browse the result in your specified output directory (an html file) using your web browser.

---
#### Start the desktop-based version 

1. On ``Linux`` or ``macOS`` systems, you may use the command line below in ``Terminal`` to open the software.

    ```
    $ cd /full/path/to/panGraphViewer/panGraphViewerApp ##Please modify the path based on your ENV
    $ python panGraphViewerApp.py
    ```
2. On ``Windows`` systems, you may search and open ``Anaconda Prompt (miniconda3)`` first and then move to the ``panGraphViewer`` directory. For example, if you have put ``panGraphViewer`` on your ``Desktop`` and the opened ``Anaconda Prompt (miniconda3)`` is in your ``C`` drive, you may use the command line below to start the program:

    ```
    > cd C:\Users\%USERNAME%\Desktop\panGraphViewer\panGraphViewerApp
    > python panGraphViewerApp.py
    ```

    If you have put ``panGraphViewer`` on another drive, you may need to move to the target drive first. For instance, the target drive is ``D``, you can move to the drive by typing **D:** in ``Anaconda Prompt (miniconda3)`` and press ``Enter`` first and then move to the ``panGraphViewerApp`` directory to execute by typing 
	```
	python panGraphViewerApp.py
	```

    Please **note** that on ``Windows`` systems, you need to use backslash ``\ `` rather than the common slash ``/`` to go to the target directory.

3. The **logging information** will show in ``Anaconda Prompt (miniconda3)`` or ``Terminal`` depending on the system that you use (this will be good for you to monitor the status of the application). 

---
### Web-based panGraphViewer

The ``web browser-based panGraphViewer`` offers administrative functions to help create accounts for different users and it allows root/super users to easily manage different accounts and enables file sharing on the Linux server. Most functions provided in the ``Desktop-based`` version have been implemented in the ``Web browser-based`` version. The ``Web-based`` application can also be installed locally if users have difficulties in installing the ``desktop-based`` application. 

#### Library installation for the web-based version
Users can use ``pip`` directly to install the needed ``python3`` libraries after moving to the ``panGraphViewerWeb`` directory.

```
pip install -r requirements.txt           ## On Linux or macOS systems
pip install -r requirements_windows.txt   ## On Windows systems
```

Users can also use ``conda`` to install the python libraries
```
conda config --add channels conda-forge
conda config --add channels bioconda
conda install django==3.1.3 django-bootstrap4 requests django-environ djangorestframework configparser natsort attrdict networkx bokeh==2.2.3 dna_features_viewer pandas whitenoise
```
As mentioned in the ``desktop-based`` version, ``pysam`` cannot be installed on ``Windows`` systems. Users need to install alternatives on ``Windows`` by using such as 
```
> conda install m2-base pyfaidx
```
For ``Linux`` or ``macOS`` users, ``pysam`` can be installed directly using 
```
$ conda install pysam
```

**Note**: ``fontawesome-free==5.15`` is needed in the ``web-based`` application, which cannot be installed using ``conda``. ``pip`` can be used to install this.

---
#### Start the web-based version
After the installation above, users can move to the ``panGraphViewerWeb`` directory by referring to the steps mentioned in the ``desktop`` version through ``Terminal`` or ``Anaconda Prompt (miniconda3)``. 

Once moving to the ``panGraphViewerWeb`` directory, users can start the application by typing 
```
python manage.py runserver <IPaddress>   ## On local machine the IPaddress can be: localhost:8000
```

or users can use the ``command`` line below to start the ``Web browser-based`` version 
```
$ bash run.sh     ## On Linux or macOS systems.
> run.bat         ## On Windows systems
```

**Note:** the IP ``0.0.0.0`` in ``run.sh`` can be modified accordingly to such as ``localhost``

Once the words ``Starting development server at http://0.0.0.0:8000/`` or similar information is showing in your ``Terminal`` or ``Anaconda Prompt (miniconda3)``, users can open a browser to open the ``web-based panGraphViewer``. For local machine, users can type ``localhost:8000`` in the browser to open panGraphViewer.

The admin page is something like ``http://localhost:8000/admin`` depending on the IP used and the initial admin info is:
```
Account:  admin
password: abcd1234
```

**Note:** please use the **go back** button provided by the web browser to **move back** rather than directly clicking the corresponding functions in the web page to perform analyses.

---

### Files accepted in the application

#### GFA files

PanGraphViewer is designed based on the reference GFA (``rGFA``) file given the flexibility of this data format. We can also accept the native [GFA_v1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) files with paths (``P``) to display graphs. However, the graph may be less informative due to the nature of GFA_v1, for example, the native GFA_v1 may lose the information of the origin of nodes.

1. If you have multiple high-quality genome assemblies from different individuals, you may use [minigraph](https://github.com/lh3/minigraph) (``Linux`` preferred) to generate an ``rGFA`` file.

    Before running, the header of the fasta file needs modifying. For example, if you have a fasta file from **Sample1** with a header like:
    ```
    >chr1
    AAAAAGCCGCGCGCGCTTGCGC
    ```

    You may modify the header to:
    ```
    >Sample1||chr1
    AAAAAGCCGCGCGCGCTTGCGC
    ```
    On ``Linux``, the command lines that can be used to achieve this are:
    ```
    $ sample=""  ## the name of the sample. For instance: Sample1
    $ fasta=""   ## full path to the fasta file
    $ name=`echo $fasta | rev | cut -d"." -f2-| rev`
    $ sed -e "s/>/>${sample}||/g" $fasta > ${name}.headerModified.fasta
    ```

    Here, we also provide a python script ``renameFastaHeader.py`` to help this conversion. The script can be found in the ``scripts`` folder under ``panGraphViewer`` --> ``panGraphViewerApp``. Or users can use the Desktop version to convert by clicking ``Tools`` --> ``Format Conversion`` --> ``Modify FASTA Header``.

    ```
    usage: renameFastaHeader.py [-h] [--version] [-f FASTA] [-n NAME] [-o OUTPUT]

    rename the header of a given fasta file

    optional arguments:
      -h, --help  show this help message and exit
      --version   show program's version number and exit
      -f FASTA    a fasta format file
      -n NAME     name of the sample
      -d DELIM    delimiter. Default: '||'
      -o OUTPUT   the output directory
    ```

    Please **note** that:

    I). If you do not modify the header of your fasta file and directly use ``minigraph`` to generate the ``rGFA`` file, ``panGraphViewer`` can still read the file, while many features, such as ``where the node comes from`` would not show in detail. A warning message will display in both UI and the opened ``Terminal`` or ``powershell``.

    II). For the sample name, please ``DO NOT`` include ``||``.

2. If you don't have an ``rGFA`` file, but a ``GFA_v1`` file with paths (``P``), you may try to follow the standard [here](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) to convert your ``GFA`` file into an ``rGFA`` file. 
   
   In case it is difficult to do so, we have provided an internal function to help convert such a ``GFA`` file to an ``rGFA`` file. Users can simply select a ``GFA`` file to browse the underlying graph, but a warning message will show if a GFA rather than an rGFA is selected. Also, our program will check internally if the input GFA_v1 has ``PATH(P)`` in it. If not, an error message will display. 
   
   Notably, if the gfa_v1 file is ``big`` (> 1Gb), our program will take a while to perform the conversion from gfa1 to rGFA. We recommend using ``gfa2rGFA.py`` located in the ``scripts`` folder to perform the conversion if the file is big.
   
---

#### VCF files

A ``VCF`` file is also accepted to show the genome graph. Basically, a reference FASTA file is **optional** if the ``VCF`` is a standard one. The program will check the input ``VCF`` file and evaluate if the ``VCF`` file meets the requirement automatically. If not, a warning or an error message will show.

Depending on your purpose, ``VCF`` **filtration** is highly recommended before plotting the underlying graph. 

Here, we also provide a method to help convert a ``VCF`` file to an ``rGFA`` file. Users can perform the conversion directly through the interface provided in the application or directly use ``vcf2rGFA.py`` under the ``panGraphViewer`` --> ``panGraphViewerApp`` --> ``scripts`` folder. 

**Note:** if there are **many** variations in the ``VCF`` file, we recommend using ``vcf2rGFA.py`` directly to convert by chromosomes/sequences rather than converting entirely. This will save a lot of computational resources when plotting graphs.

The usage of ``vcf2rGFA.py`` is shown below. Both ``Windows`` and ``Linux/macOS`` users can directly use this script to convert a ``VCF`` file to an ``rGFA`` file. 
```
usage: vcf2rGFA.py [-h] [--version] [-f FASTA] [-b BACKBONE] [-v VCF] [-o OUTPUT] [-c [CHR [CHR ...]]] [-n NTHREAD]
    
Convert a vcf file to an rGFA file
    
optional arguments:
    -h, --help          show this help message and exit
    --version           show program's version number and exit
    -f FASTA            a fasta format file that from the backbone sample
    -b BACKBONE         the name of the backbone sample
    -v VCF              the vcf file
    -o OUTPUT           the output directory
    -c [CHR [CHR ...]]  the name of the chromosome(s) [default: all chroms]
    -n NTHREAD          number of threads [default: 4]
```
---

#### Annotation files 

If users want to check nodes that fall in a gene model region, a ``BED``, ``GTF`` or ``GFF3`` file from the ``backbone`` sample can be provided to do so. Basically, the ``BED`` file should contain at least 6 columns as shown below.

| Column  | Information                                                           |
| :-----: | :---------------------------------------------------------------------|
| 1       | Chromosome ID                                                         |
| 2       | Gene start position                                                   |
| 3       | Gene end position                                                     |
| 4       | Gene ID                                                               |
| 5       | Score (or others; the program does not use the info in this column)   |
| 6       | Orientation                                                           |

* Users can load the ``BED``, ``GTF`` or ``GFF3`` file to the application to check the coordinate overlaps between nodes and genes. 
* By default, genes overlapping with more than ``2`` nodes will be shown in the dropdown menu. 
* A gene list will be saved in the output directory after parsing the annotation file.
* When using this function, a plot containing the selected gene and nodes falling in the gene region will be shown. A subgraph of the gene region will also be shown.

---

### Q&A:

#### Operating system 

* For the ``desktop-based`` application, it has been optimized on 
	```
	Windows 10
	macOS Big Sur, macOS Monterey and
	Ubuntu 18.04.5 
	```
	
	For other operating systems or equivalents, our tool may also work. However, on older operating systems, such as ``Ubuntu 16.04``, ``PyQtWebEngine`` may not work properly. 

* For the ``web browser-based`` version, we suggest running in ``Linux`` or ``macOS`` environment. If users want to run on ``Windows`` systems, ``Windows 10`` or above is recommended. 

	Users can also use [docker](https://docs.docker.com/get-docker/) to run the ``web browser-based`` version (see the instruction [here](panGraphViewerWeb/README.md)). However, ``WSL`` is needed to run the docker version on ``Windows 10`` or above.

---
#### Minimum computational resource needed
```
Memory:  1Gb
Threads: 2
```
Users can adjust the RAM depending on the size of the ``rGFA`` file. Sometimes, the RAM used at some point could be equal to or slightly above the size of the ``rGFA`` file.

---
#### Which application should I use 

* Depending on the purpose and preference, users can select either the ``desktop-based`` application or the ``web-based`` application. Basically, the ``desktop-based`` application is more intendedly designed for single users who have little computer science background. This application allows users to easily browse the graph by simply clicking the buttons. 

* By contrast, the ``web-based`` application is more suitable for multiple users who have a shared Linux server/cluster. The administrator can easily deploy the application on the server and set accounts for different users to use this tool. Users can also easily share their files or results on the server. While the ``web-based`` application is not limited to a server, it can be installed locally in case that the ``desktop-based`` version fails in installation.

* Most of the functions in the ``desktop-based`` version and the ``web-based`` version are the same. However, there are some differences. First, to make the response time acceptable, in the ``web-based`` application, we only use ``cytoscape.js`` to draw graphs. In contrast, ``vis.js`` and ``cytoscape.js`` are used in the ``desktop-based`` application to draw graphs depending on how many nodes a user wants to browse. By default, ``200`` nodes are allowed to show in ``vis.js``-based graph. Users can change this setting in ``Settings``--> ``Graph Modification`` to decide either using ``vis.js`` or ``cytoscape.js``. The reason to use two different javascript libraries is that ``vis.js`` can draw a nicer graph, but needs more loading time if there are tens of thousands of nodes. While ``cytoscape.js`` does not have such a problem. Moreover, in the ``web-based`` application, users can easily browse the graph and node sequences without the problem of system differences.

---
#### How to use the application

* For more detailed steps to run ``panGrapViewer``, please refer to the [Manual](doc/Manual.md)

---
#### Backbone sample

* The ``backbone`` sample is the one used as the **main sequence provider** to produce the pangenome graph or the reference sample to produce the ``VCF`` file. 
* In a pangenome graph, most of the nodes are from the ``backbone`` sample (shared by all) with some nodes (variations) from other samples. 

---
#### Colors shown in the graph

* Each sample uses one particular colour and **the most frequent colour** is the one used for the ``backbone`` sample. 
* By default, the colours are randomly selected by the program from a designed colour palette. 
* Users can customise the color scheme in the ``config_default.ini`` or ``config.ini`` file to fix the colors for samples (number of colors = number of samples).
* In the ``legend``, users can check the color used for each sample.

---
#### Type of graphs 

* We provide two kinds of graph plots in the ``desktop-based`` application to achieve a good performance and visualisation. By default, if the number of checked nodes <= **200**, ``vis.js`` based graph will show. Otherwise, a ``cytoscape.js`` based graph will show. Users can change the settings in the ``desktop-based`` application to set the type of plot.
* In the ``web-based`` application, ``cytoscape.js`` based graph is provided.

---
#### Shapes shown in the graph

* If users use a ``VCF`` file to show graphs, we use different nodes shapes to represent different kinds of variants. For instance, in the default settings for the ``vis.js``-based graph,
	```
	o represents SNP
	△ represents deletion
	▽ represents insertion 
    ☷ represents duplication 
	text shows inversion
	☆ represents translocation
	```

Users can change the corresponding settings to select preferred node shapes to represent different variations on the ``desktop-based`` application.

In the ``legend``, users can check the shapes used for representing different variations.

---
#### Different variations 

* If users use a ``VCF`` file to generate a graph genome, when moving the mouse to the graph node, the program will show the variation types automatically, such as 
	```
	SNP:   single nucleotide polymorphism
	INS:   insertion
    DEL:   deletion
	INV:   inversion 
	DUP:   duplication
    TRANS: translocation
	``` 

The corresponding nodes from the ``backbone`` sample will also be linked and shown. 

---
#### Pros and cons of the application

There are some pros and cons of this application. We list some here for your reference 

***Pros***

* can directly plot a subgraph of interest from a pangenome graph
* can specify coordinates to check subgraphs and nodes 
* can browse graphs from a VCF file 
* can show variations, particularly structural variations if the graph is generated from a VCF file
* can check nodes falling in a gene model region and alter users if the nodes have the chance to change the function of the gene model in some individuals   
* can show the origin of particular nodes in a pangenome graph
* Computational resource-efficient in using this tool

***Cons***

* When checking hundreds or thousands of nodes in vis.js-based graph, the nodes may cluster together in the plot which may be difficult to interpret [Third-party javascript library layout problem; in this case we suggest switching to cytoscape.js-based graphs in Settings]
* Some python libraries used, such as 'PyQtWebEngine' may not be compatible with all operating systems if using the desktop-based application [Third-party python library problem; in this case users may use the web browser-based version]
* Remote accessing to the desktop-based application may confront plotting problems, particularly when using macOS system remote accessing to the Linux system [Selection problem; in this case users may want to install the desktop-based application locally or use the web browser-based version]
* In macOS, if users use the web browser-based application, the hover for nodes cannot show properly in ``Safari`` [Browser problem; in this case users may use ``Chrome`` or ``Microsoft Edge``]
* May not be able to display a graph containing more than ``20,000`` nodes in the Desktop-based application and ``5,000`` nodes in the Web browser-based application [Rendering problem - may render slowly or even crash in the display canvas. This is a problem inherited from the third-party javascript libraries we used. In this case, users may want to adjust the plotting region to avoid this problem]

---
Enjoy using panGraphViewer!
