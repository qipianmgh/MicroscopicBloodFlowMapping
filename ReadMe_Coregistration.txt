Instruction on TPM & OCT data co-registration v1.0 
1. Prepare angiogram files of TPM and OCT. 
   TPM_Stack_filtered_05022018.mat
   OCT_DLS_angio_05022018.mat
2. Run coregisterOCTand2PM.m from file path for 3D global co-registration. 
3. Generate mini-stacks of 40 um depth, then run the coregisterOCTand2PM.m again for each mini-stack pairs of TPM and OCT to perform 2D co-registration. 
4. Retrieve the transformation matrices from both 3D and 2D co-registration.
5. Read TPM graphs for re-graph and manual correction in VIDA to process the graphs. 
6. Perform co-registration on the TPM graphs using the transformation matrices from Step 4.
7. Feed the results into ‘Mapping_OCT_Flow.m’ for flow speed mapping. 
