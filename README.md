Overview
--------

This toolbox provides support utilities for the Photorealistic Virtual City dataset.

http://people.csail.mit.edu/biliana/projects/iccv2011/

We provide additional tools that define image to image correspondences for all pairs
of overlapping images.


Installation
------------

This toolbox relies on the nested shape descriptor toolbox for displaying correspondences and the original PVC dataset. 

1. Download [PVC dataset](http://people.csail.mit.edu/biliana/projects/iccv2011/) and unpack into $indir

2. Download and install nested shape descriptor toolbox

sh\> git clone https://github.com/jebyrne/seedoflife  
\>\> cd seedoflife    
\>\> set_paths    


3. Download [precomputed PVC correspondences](http://dropbox.com) and unpack to ${outdir}

These correspondences can also be generated manually:

\>\> cd ${PhotorealisticVirtualWord}    
\>\> pvc_correspondence(${indir}, ${outdir}    

4. Run demo    
\>\> cd ${PhotorealisticVirtualWord}    
\>\> demo_pvc(${indir}, ${outdir}    

5. To run a descriptor comparison and plot performance results 

\>\> cd ${PhotorealisticVirtualWord}    
\>\> eval_pvc(${indir}, ${outdir}    


References
----------
* Biliana Kaneva, Antonio Torralba, William Freeman. "[Evaluation of Image Features Using a Photorealistic Virtual World](http://people.csail.mit.edu/biliana/projects/iccv2011/). ICCV, 2011
* J. Byrne and J. Shi, "[Nested Shape Descriptors](https://www.dropbox.com/s/g8yc76ffx8ia99d/iccv13_nsd_final.pdf)", International Conference on Computer Vision (ICCV'13), Sydney Australia, 2013
* Engin Tola, Vincent Lepetit, Pascal Fua, [DAISY](http://cvlab.epfl.ch/software/daisy): An Efficient Dense Descriptor Applied to Wide Baseline Stereo, IEEE Transactions on Pattern Analysis and Machine Intelligence Vol. 32, Nr. 5, pp. 815 - 830, May 2010


