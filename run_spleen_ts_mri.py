import pandas as pd
import sys, glob, os, fnmatch, shutil, subprocess, csv, pydicom
import SimpleITK as sitk
from scipy import misc
import scipy
import math
import matplotlib.image as mp
from shutil import copyfile
import pandas as pd
import numpy as np
import shutil
import collections
import re
import nibabel as nib
from totalsegmentator.python_api import totalsegmentator
from totalsegmentator.map_to_binary import class_map
from totalsegmentator.dicom_io import *

import dicom2nifti



os.environ["CUDA_VISIBLE_DEVICES"] = "1"


def sort_nicely(l):
    l.sort(key=alphanum_key)
    
    
def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    return [tryint(c) for c in re.split('([0-9]+)', s) ]

















if __name__ == "__main__":

    #Input directory for drive
    input_dir = "/mypool/shared/MRI/MGBB"
        
    #Output directory for segmentations
    output_dir = "/mypool/shared/MRI/MGBB_Segmentations/"


    data_path = "Spleen_Unique_Series_Only.csv"

    df = pd.read_csv(data_path)
     
    #Organs to keep
    organs_from_full = ['spleen']

    #First iteration
    first = True

    nifti_dir = "tmp_nifti"
    if not os.path.exists(nifti_dir):
        os.mkdir(nifti_dir)

    #Clear tmp nifti directory
    nft_files = os.listdir(nifti_dir)
    for f in nft_files:
        os.remove(os.path.join(nifti_dir,f))

    all_studies = pd.unique(df.ORIG_SUID)

    #Go through the 165 unique study descriptions -> hopefully patterns in the sequences obtained
    #Come up with a systematic way to analyze each of these studies (including stitching, identifying in/opp, running TS, running marco's)

    #TODO - make the loop per study instead of row
    #Handle cases where in/opp are split
    #Only run marco bodycomp on dixon
    #Run TotalSeg on everything
    #Lastly some cases have no mention of in/opp (dualecho)
    #And some cases have fat/water series as well - need to handle to calculate SMFF
    #Finally can we generate fat/water using in/opp and acquisition parameters

    spleen = np.zeros(df.shape[0])
    sat = np.zeros(df.shape[0])
    vat = np.zeros(df.shape[0])
    muscle = np.zeros(df.shape[0])
    completed = np.zeros(df.shape[0])
    errors = []
    for i in range(df.shape[0]):
        print("MR: " + str(i) + ", out of " + str(df.shape[0]))
        path = os.path.join(input_dir,df['Path'][i])
        key = df['Path'][i].replace("/","_")
        fname = key + "_Segmented.nii"
        #Get unique echo time values and separate these into separate folders
        dicoms = os.listdir(path)
        
        
        all_dicoms = [pydicom.read_file(os.path.join(path,f)) for f in dicoms]
        try:
            echo_times = [float(f[0x18,0x81].value) for f in all_dicoms]
            series_numbers = [int(f[0x20,0x11].value) for f in all_dicoms]
            acqs = [int(f[0x20,0x12].value) for f in all_dicoms]
            inst = [int(f[0x20,0x13].value) for f in all_dicoms]
        except:
            errors.append("Missing DICOM tag")
            print("Missing DICOM tags")
            continue

        
        
        #Get files corresponding to current te
        
        #Move files to tmp_dicom directory
        
        #convert that directory
        

        #clear directory
        if len(np.unique(echo_times))==2: #In and opposed-phase
            tmp_dicom = "tmp_dicom"
            if not os.path.exists(tmp_dicom):
                os.mkdir(tmp_dicom)
            else:
                shutil.rmtree(tmp_dicom)
                os.mkdir(tmp_dicom)
            
            unique_times = np.sort(np.unique(echo_times))
            
            niftis = []
            try:
                for t in unique_times:
                    indices  = np.where(echo_times==t)[0]
                    curr_files = [dicoms[i] for i in indices]
                    for cf in curr_files:
                        shutil.copy(os.path.join(path,cf),os.path.join(tmp_dicom,cf))

                    dicom2nifti.convert_directory(tmp_dicom,nifti_dir)
                    out = nib.load(os.path.join(nifti_dir,os.listdir(nifti_dir)[0]))
                    niftis.append(out)
                    shutil.rmtree(tmp_dicom)
                    os.mkdir(tmp_dicom)
            except:
                errors.append("NIFTI conversion in/opp")
                print("Error converting in/opp to NIFTI")
                continue
        #Niftis[0] is in-phase
        #Niftis[1] is opposed-phase
        
            out = niftis[0]
        
        else:      
        ###Convert dicom to nifti as done on RAP
            try:
                dicom2nifti.convert_directory(path, nifti_dir)
                out = nib.load(os.path.join(nifti_dir,os.listdir(nifti_dir)[0]))
            except:
                print("Error converting DICOM to NIFTI")
                errors.append("NIFTI conversion normal")
                continue


        #Add logic to decide whether to stitch, whether to directly load nifti, whether to do any other adjustments

        #Figure out the format of "out" and save to file along with a nifti of the original MRI
        #import pdb; pdb.set_trace()
        if(not os.path.exists(os.path.join(output_dir,fname))):
            try:
                #Running out of memory during totalsegmentator call -- why???
                #Then, run total_mr and save the spleen
                segmented = totalsegmentator(out,task="total_mr")
                vox_volume = np.prod(segmented.header.get_zooms())
                spleen[i] = len(np.where(segmented.get_fdata()==1)[0])*vox_volume
                nib.save(segmented,os.path.join(output_dir,fname))
            except:
                errors.append("TS error")
                print("Issue with total segmentation for: " + str(path))
                spleen[i] = -1
        else:
            segmented = nib.load(os.path.join(output_dir,fname))
            vox_volume = np.prod(segmented.header.get_zooms())
            spleen[i] = len(np.where(segmented.get_fdata()==1)[0])*vox_volume
        if(not os.path.exists(os.path.join(output_dir,key+"_TS_BC.nii"))):
            try:
                seg_tissue = totalsegmentator(out,task="tissue_types_mr")
                sat[i] = np.prod(seg_tissue.header.get_zooms())*len(np.where(seg_tissue.get_fdata()==1)[0])
                vat[i] = np.prod(seg_tissue.header.get_zooms())*len(np.where(seg_tissue.get_fdata()==2)[0])
                muscle[i] = np.prod(seg_tissue.header.get_zooms())*len(np.where(seg_tissue.get_fdata()==3)[0])
                #Save full segmentation mask (out)
                nib.save(seg_tissue,os.path.join(output_dir,key+"_TS_BC.nii"))
            except:
                errors.append("TS Body Comp Error")
                print("Issue with fat/muscle segmentation for: " + str(path))
                sat[i] = -1
                vat[i] = -1
                muscle[i] = -1
        else:
            seg_tissue = nib.load(os.path.join(output_dir,key+"_TS_BC.nii"))
            sat[i] = np.prod(seg_tissue.header.get_zooms())*len(np.where(seg_tissue.get_fdata()==1)[0])
            vat[i] = np.prod(seg_tissue.header.get_zooms())*len(np.where(seg_tissue.get_fdata()==2)[0])
            muscle[i] = np.prod(seg_tissue.header.get_zooms())*len(np.where(seg_tissue.get_fdata()==3)[0])
        print("Completed segmentation")

        #Grab x,y,z dimensions

        
        print("Spleen: " + str(spleen[i])) 
        print("SAT: " + str(sat[i])) 
        print("VAT: " + str(vat[i])) 
        print("Muscle: " + str(muscle[i]))
        completed[i] = 1
        errors.append("None")
        
        
        
 
        if(not os.path.exists(os.path.join(output_dir,key+"_Original.nii"))):
            #Save original MRI image cube
            try:
                nib.save(out,os.path.join(output_dir,key+"_Original.nii"))
            except:
                errors[len(errors)-1]  = "Issue saving nifti"
        
        nft_files = os.listdir(nifti_dir)
        for f in nft_files:
            os.remove(os.path.join(nifti_dir,f))
        print("Done")
    df['Splenic_Volume'] = spleen
    df['SAT'] = sat
    df['VAT'] = vat
    df['Muscle'] = muscle
    df['Error'] = errors
    df.to_csv(data_path.replace(".csv","_Segmented.csv"),index=False)
