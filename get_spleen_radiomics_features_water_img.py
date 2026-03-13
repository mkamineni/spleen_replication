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
from radiomics import featureextractor
from fastprogress import progress_bar
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




def getRadiomics(nib_image,mask_image):
    image_data = nib_image.get_fdata()
    mask_data = mask_image.get_fdata()
    # Step 3: Get the image metadata/affine information
    affine = nib_image.affine
    header = nib_image.header
    pixdim = header.get_zooms()  # Voxel dimensions

    # Step 4: Convert numpy arrays to SimpleITK images (required by pyradiomics)
    sitk_image = sitk.GetImageFromArray(image_data)
    sitk_mask = sitk.GetImageFromArray(mask_data)

    # Step 5: Set the physical spacing information (important for correct feature calculation)
    spacing = (float(pixdim[0]), float(pixdim[1]), float(pixdim[2]))

    # Set the spacing properly
    sitk_image.SetSpacing(spacing)
    sitk_mask.SetSpacing(spacing)

    # Step 6: Initialize the feature extractor
    extractor = featureextractor.RadiomicsFeatureExtractor()
    # Optional: Customize settings
    # extractor.settings.update({'binWidth': 25, 'resampledPixelSpacing': [1, 1, 1]})
    
    

    # Disable all features first
    extractor.disableAllFeatures()

    # Then enable only the required feature classes
    extractor.enableFeatureClassByName('firstorder')
    extractor.enableFeatureClassByName('shape')
    extractor.enableFeatureClassByName('glcm')
    extractor.enableFeatureClassByName('glrlm')
    extractor.enableFeatureClassByName('glszm')
    extractor.enableFeatureClassByName('gldm')

    # Disable all features within these classes
    extractor.enableFeaturesByName(firstorder=[], shape=[], glcm=[], glrlm=[], glszm=[], gldm=[])

    # Then enable only the specific features you want
    # The feature names need to match pyradiomics' naming convention
#    extractor.enableFeaturesByName(
#        firstorder=['Energy'],
#        shape=['Sphericity'],
#        glcm=['Id','Idn', 'Idmn', 'Correlation'],
#        glrlm=['RunLengthNonUniformity'],
#        glszm=['LargeAreaLowGrayLevelEmphasis', 'GrayLevelNonUniformity'],
#        gldm=['SmallDependenceHighGrayLevelEmphasis', 'GrayLevelVariance']
#    )

    extractor.enableFeaturesByName(glcm=['Id'])
    # Add these checks before extraction:
    print(f"Image size: {sitk_image.GetSize()}")
    print(f"Mask size: {sitk_mask.GetSize()}")
    print(f"Image spacing: {sitk_image.GetSpacing()}")
    print(f"Mask spacing: {sitk_mask.GetSpacing()}")

    # Check if mask has valid regions
    mask_array = sitk.GetArrayFromImage(sitk_mask)
    print(f"Mask value range: {mask_array.min()} to {mask_array.max()}")
    print(f"Non-zero voxels: {np.count_nonzero(mask_array)}")
    import pdb; pdb.set_trace()
    try:
        featureVector = extractor.execute(sitk_image, sitk_mask)
    except:
        return None
    # Extract features
    #featureVector = extractor.execute(sitk_image, sitk_mask)

    # Create a list to store the features and values
    feature_list = []

    # Map from pyradiomics feature names to your desired display names
    feature_name_map = {
        'original_firstorder_Energy': 'Energy',
        'original_shape_Sphericity': 'Sphericity',
        'original_glcm_Idm' : 'IDM',
        'original_glcm_Idn' : 'IDN',
        'original_glcm_Idmn' : 'IDMN',
        'original_glcm_Id' : 'ID',
        'original_glcm_InverseDifference': 'GLCM inverse difference',
        'original_glcm_Correlation': 'GLCM correlation',
        'original_glrlm_RunLengthNonUniformity': 'GLRLM run length non-uniformity',
        'original_glszm_LargeAreaLowGrayLevelEmphasis': 'GLSZM large area low gray-level emphasis',
        'original_glszm_GrayLevelNonUniformity': 'GLSZM gray-level non-uniformity',
        'original_gldm_SmallDependenceHighGrayLevelEmphasis': 'GLDM small dependence high gray-level emphasis',
        'original_gldm_GrayLevelVariance': 'GLDM gray-level variance'
    }

    # Create organized output
    results = {}
    print('\nExtracted features:')
    for featureName, featureValue in featureVector.items():
        # Only process the features we're interested in
        if featureName in feature_name_map:
            display_name = feature_name_map[featureName]
            results[display_name] = float(featureValue)
            print(f"{display_name}: {featureValue}")
            feature_list.append([display_name, featureValue])

    return(results)



















if __name__ == "__main__":

    #Input directory for DICOM files
    input_dir = "/mypool/shared/MRI/MGBB"

    #Output directory for segmentations
    output_dir = "/mypool/shared/MRI/MGBB_Segmentations/"

    #.csv with list of files to parse
    data_path = "Re_Run_No_In_And_Opposed_072325.csv"
    output_data_path = "Second_Round_No_Water_ID_082425.csv"
   

    nifti_dirs = ["opp_phase","in_phase"]
    

    df = pd.read_csv(data_path)
     
    #First iteration
    first = True
    fnames = []
    
    for i in range(df.shape[0]):
        #import pdb; pdb.set_trace()
        #if (df.Dixon_Water_Radiomics[i]==1) or (df.Dixon_Water_Radiomics[i]==2):
        #    continue
        print("MR: " + str(i) + ", out of " + str(df.shape[0]))
        path = os.path.join(input_dir,df['Path'][i])
        key = df['Path'][i].replace("/","_")
        fname = key + "_Segmented.nii"
        #Get unique echo time values and separate these into separate folders
        dicoms = os.listdir(path)
        
        for nd in nifti_dirs:
            if not os.path.exists(nd):
                os.mkdir(nd)
            nft_files = os.listdir(nd)
            for f in nft_files:
                os.remove(os.path.join(nd,f))
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
            
            
            #Opposed phase is shorter echo time
            niftis = []
            try:
                for t in range(len(unique_times)):
                    indices  = np.where(echo_times==unique_times[t])[0]
                    curr_files = [dicoms[i] for i in indices]
                    for cf in curr_files:
                        shutil.copy(os.path.join(path,cf),os.path.join(tmp_dicom,cf))

                    dicom2nifti.convert_directory(tmp_dicom,nifti_dirs[t])
                    out = nib.load(os.path.join(nifti_dirs[t],os.listdir(nifti_dirs[t])[0]))

                    niftis.append(out)
                    shutil.rmtree(tmp_dicom)
                    os.mkdir(tmp_dicom)
                
                opp_p = niftis[0]
                in_p = niftis[1]
                
                water = (in_p.get_fdata() + opp_p.get_fdata())/2
                
                water_nii = nib.Nifti1Image(water,in_p.affine)

                
            except:
                df.loc[i,"Dixon_Water_Radiomics"] = -1
                print("Error converting in/opp to NIFTI")
                continue
        #Niftis[0] is opposed-phase
        #Niftis[1] is in-phase
        elif len(np.unique(echo_times))>2:
            df.loc[i,"Dixon_Water_Radiomics"] = -2 # more than 2 echo times
        else:
            df.loc[i,"Single_Echo_Time"] = echo_times[0]
            df.loc[i,"Dixon_Water_Radiomics"] = 2 #Not water phase but something else - single echo time
            continue

        import pdb; pdb.set_trace()

        #Assuming segmentation was successful since we have a splenic volume in range
        if(os.path.exists(os.path.join(output_dir,df.Filename[i]))):
            segmented = nib.load(os.path.join(output_dir,df.Filename[i]))
            
            # Convert to binary (should already be binary but to ensure)
            spleen_data = segmented.get_fdata()
            
            #convert to binary
            spleen_data[spleen_data!=1] = 0
            
            # Create a new binary mask
            binary_spleen = nib.Nifti1Image(spleen_data.astype(np.int8), segmented.affine)
            features = getRadiomics(water_nii,binary_spleen)
            
            # Update the DataFrame row with features
            if features is not None:
                df.loc[i, features.keys()] = features.values()
                df.loc[i,"Dixon_Water_Radiomics"] = 1
            
        if i%1000==0:
            df.to_csv(output_data_path.replace(".csv","_" + str(i) + ".csv"),index=False)
    df.to_csv(output_data_path,index=False)    
