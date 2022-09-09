# PanGraphViewer -- show panenome graphs in an easy way

![PyPI - Python Version](https://img.shields.io/pypi/pyversions/Django?color=green&logoColor=green&style=social) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/tf-chan-lab/panGraphViewer?color=green&logoColor=green&style=social) ![GitHub](https://img.shields.io/github/license/tf-chan-lab/panGraphViewer?color=green&logoColor=green&style=social) ![GitHub all releases](https://img.shields.io/github/downloads/tf-chan-lab/panGraphViewer/total?color=green&logoColor=green&style=social) 

### Table of Contents
<img align="right" width="300" src="doc/src/demo.jpg"/>

- [Versions and installation](#versions-and-installation)
- [Files accepted ](#files-accepted)
- [Start the program](#start-the-program)
- [Run the program](#run-the-program)
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
## Versions and installation
Here we provide **two** application versions:
```
● Desktop-based application
● Web-based application
```
Overall, **Python3.6 or above** is needed to run this software. We highly recommend users using ``miniconda3`` to run this program.

Users can refer to the [Installation](https://github.com/TF-Chan-Lab/panGraphViewer/wiki/Installation) at the [Wiki](https://github.com/TF-Chan-Lab/panGraphViewer/wiki) page to install the tool.

To futher assist the installation and use, we have provided two ``Youtube videos`` to guild users to use this tool.

**On Windows** 

Users may refer to https://youtu.be/YrltzD8R5Io

**On Other operating systems**

Users may refer to https://youtu.be/tpLjcOz2E4U

## Files accepted 

``PanGraphViewer`` accepts different pangenome graph formats, including ``rGFA``, ``GFA_v1`` and ``VCF``. 

``PanGraphViewer`` can also accept genome annotation files, such as ``BED``, ``GTF/GFF`` files.

Before using this tool, users may refer to [data formats](https://github.com/TF-Chan-Lab/panGraphViewer/wiki/Files-accepted) to prepare the files and then run the program.

## Start the program
Users may refer to [start the program](https://github.com/TF-Chan-Lab/panGraphViewer/wiki/Start-the-program) to open the user interface.

## Run the program 
Users may refer to the [Manual](https://github.com/TF-Chan-Lab/panGraphViewer/blob/main/doc/Manual.md) to run the program.

## Q&A:

### Operating system 

* For the ``desktop-based`` application, it has been optimized on 
	```
	Windows 10
	macOS Big Sur, macOS Monterey, Ubuntu 18.04.5 and Ubuntu 20.04
	```
	
	For other operating systems or equivalents, the tool may also work. However, on older operating systems, such as ``Ubuntu 16.04``, ``PyQtWebEngine`` may not work properly. 

* For the ``web-based`` version, we suggest running in ``Linux`` or ``macOS`` environment. If users want to run on ``Windows`` systems, ``Windows 10`` or above is recommended. 

	Users can also use [docker](https://docs.docker.com/get-docker/) to run the ``web-based`` version (see the instruction [here](panGraphViewerWeb/README.md)). However, ``WSL`` is needed to run the docker version on ``Windows 10`` or above.

---
### Minimum computational resource needed
```
Memory:  2Gb
Threads: 2
```
Users can adjust the RAM depending on the size of the pangenome graph. 

From our preliminary tests, it seems ``8Gb`` RAM would be sufficient to view most graphs.

---
### Which application should I use 

* Depending on the purpose and preference, users can select either the ``desktop-based`` application or the ``web-based`` application. Basically, the ``desktop-based`` application is more intendedly designed for single users who have little computer science background. This application allows users to easily browse the graph by simply clicking the buttons. 

* By contrast, the ``web-based`` application is more suitable for multiple users who have a shared Linux server/cluster. The administrator can easily deploy the application on a server and set accounts for different users to  access this tool in any device with a web browser and internet connection. Users can also easily share their files or results on the server. While the ``web-based`` application is not limited to a server, it can be installed locally in case that the ``desktop-based`` version fails in installation.

* Most of the functions in the ``desktop-based`` version and the ``web-based`` version are the same. However, there are some differences. First, to make the response time acceptable, in the ``web-based`` application, we only use ``cytoscape.js`` to draw graphs. In contrast, ``vis.js`` and ``cytoscape.js`` are used in the ``desktop-based`` application to draw graphs depending on how many nodes a user wants to browse. By default, ``200`` nodes are allowed to show in ``vis.js``-based graph. Users can change this setting in ``Settings``--> ``Graph Modification`` to decide either using ``vis.js`` or ``cytoscape.js``.

---
### Backbone sample

* The ``backbone`` sample is the one used as the **main sequence provider** to produce the pangenome graph or the reference sample to produce the ``VCF`` file. 
* In a pangenome graph, most of the nodes are from the ``backbone`` sample (shared by all) with some nodes (variations) from other samples. 

---
### Colors shown in the graph

* Each sample uses one particular colour and **the most frequent colour** is the one used for the ``backbone`` sample. 
* By default, the colours are randomly selected by the program from a designed colour palette. 
* Users can customise the color scheme in the ``config_default.ini`` or ``config.ini`` file to fix the colors for samples (number of colors = number of samples).
* In the ``legend``, users can check the color used for each sample.

---
### Type of graphs 

* We provide two kinds of graph plots in the ``desktop-based`` application to achieve a good performance and visualisation. By default, if the number of checked nodes <= **200**, ``vis.js`` based graph will show. Otherwise, a ``cytoscape.js`` based graph will show. Users can change the settings in the ``desktop-based`` application to set the type of plot.
* In the ``web-based`` application, ``cytoscape.js`` based graph is provided.

---
### Shapes shown in the graph

* If users use a ``VCF`` file to show variation-based graphs, we use different nodes shapes to represent different kinds of variants. For instance, in the default settings for the ``vis.js``-based graph,
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
### Pros and cons of the application

There are some pros and cons of this application. We list some here for your reference 

***Pros***

* can directly plot a subgraph of interest from a pangenome graph
* can specify coordinates to check subgraphs and nodes 
* can browse graphs from a ``VCF`` file 
* can show variations, particularly structural variations if the graph is generated from a ``VCF`` file
* can check nodes falling in a gene model region and alter users if the nodes have the chance to change the function of the gene model in some individuals   
* can show the origin of particular nodes in a pangenome graph
* Computational resource-efficient in using this tool

***Cons***

* When checking hundreds or thousands of nodes in vis.js-based graph, the nodes may cluster together in the plot which may be difficult to interpret [Third-party javascript layout problem; in this case we suggest switching to cytoscape.js-based graphs in Settings]
* Some python libraries used, such as 'PyQtWebEngine' may not be compatible with all operating systems if using the desktop-based application [Third-party library problem; in this case users may use the web browser-based version]
* Remote accessing to the desktop-based application may confront plotting problems, particularly when using macOS system remote accessing to the Linux system [Selection problem; in this case users may want to install the desktop-based application locally or use the web browser-based version]
* In macOS, if users use the web browser-based application, the hover for nodes cannot show properly in ``Safari`` [Browser problem; in this case users may use ``Chrome`` or ``Microsoft Edge``]
* May not be able to efficiently display a graph containing more than ``20,000`` nodes [Render problem - may render slowly or even crash in the display canvas. This is a problem inherited from the third-party javascript we used. In this case, users may want to adjust the plotting region to avoid this problem]

---
Enjoy using panGraphViewer!
