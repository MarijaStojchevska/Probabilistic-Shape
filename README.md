# Reconstruction of Missing Parts of Given Partial Bones 
<p><div align="justify">
In this project, we consider 3D meshes and statistically evaluate them for further predictions using the Software Library Scalismo developed by the Graphics and Vision Research Group at the University of Basel. More precisely, we create a statistical shape model for healthy femur bones (the longest bones in the human body) using an appropriate training set. By creating the model we generate a good representation of healthy bones and use this information to complete the missing parts of the samples with missing pieces.
</div></p>

<p><div align="justify">
We perform this task stepwise: first, we apply rigid alignment to normalize the 50 training instances available on the SICAS Medical Image Repository, then we create a Gaussian Process model with smooth shape deformations for establishing correspondence among all the data sets with the help of the Iterative Closest Point method, and finally, we reconstruct the missing parts of the partial bones by creating a PCA model of the bone using the training data in correspondence.
</div></p>

More details on this project can be found in the attached report.
