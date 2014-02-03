Overview
--------

This toolbox provides support utilities for the Photorealistic Virtual City dataset.

http://people.csail.mit.edu/biliana/projects/iccv2011/

We provide additional tools that define image to image correspondences for all pairs
of overlapping images.

Correspondences between image pairs are stored in matlab MAT files, following the PVC filename convention

    asgn_{CameraIndex}_{Location2Location}_{orientation}.mat
    asgn_{CameraIndex}_{LocationIndex}_{orientation2orientation}.mat
    asgn_{CameraIndex}_{Location2Location}_{orientation2orientation}.mat

Pixel correspondences between a reference image (im_ref) and an observed image (im_obs) are represented as frames such that the kth pixel (ij) stored in matrix fr_obs(1:2,k) in image im_obs corresponds to kth pixel (ij) stored in mat.fr_obs2ref(1:2,k) in image im_ref. 

Installation
------------

This toolbox relies on the nested shape descriptor toolbox for displaying correspondences and the original PVC dataset. 

1. Download [PVC dataset](http://people.csail.mit.edu/biliana/projects/iccv2011/) and unpack into $indir

2. Download and install nested shape descriptor toolbox

        sh> git clone https://github.com/jebyrne/seedoflife  
        >> cd seedoflife    
        >> set_paths    

3. Download [precomputed PVC correspondences](https://www.dropbox.com/s/va8x0qdgfpta2xb/PhotorealisticVirtualCity.zip) and unpack to ${outdir}.  Alternatively, these correspondences can also be generated manually without downloading the ZIP file.

        >> cd ${PhotorealisticVirtualWord}    
        >> pvc_correspondence(${indir}, ${outdir})

4.  Run demo to show corresponding pixels.

        >> cd ${PhotorealisticVirtualWord}    
        >> demo_pvc(${indir}, ${outdir})

5. To run a descriptor comparison and plot performance results 

        >> cd ${PhotorealisticVirtualWord}    
        >> eval_pvc(${indir}, ${outdir})    


References
----------
* Biliana Kaneva, Antonio Torralba, William Freeman. "[Evaluation of Image Features Using a Photorealistic Virtual World](http://people.csail.mit.edu/biliana/projects/iccv2011/). ICCV, 2011
* J. Byrne and J. Shi, "[Nested Shape Descriptors](https://www.dropbox.com/s/g8yc76ffx8ia99d/iccv13_nsd_final.pdf)", International Conference on Computer Vision (ICCV'13), Sydney Australia, 2013
* Engin Tola, Vincent Lepetit, Pascal Fua, [DAISY](http://cvlab.epfl.ch/software/daisy): An Efficient Dense Descriptor Applied to Wide Baseline Stereo, IEEE Transactions on Pattern Analysis and Machine Intelligence Vol. 32, Nr. 5, pp. 815 - 830, May 2010


