# Reconstruction of Missing Parts of Given Partial Bones 
<p><div align="justify">
In this project, we consider 3D meshes and statistically evaluate them for further predictions using the Software Library Scalismo developed by the Graphics and Vision Research Group at the University of Basel. More precisely, we create a statistical shape model for healthy femur bones (the longest bones in the human body) using an appropriate training set. By creating the model we generate a good representation of healthy bones and use this information to complete the missing parts of the samples with missing pieces.
</div></p>

<p><div align="justify">
We perform this task stepwise: first, we apply rigid alignment to normalize the 50 training instances available on the SICAS Medical Image Repository, then we create a Gaussian Process model with smooth shape deformations for establishing correspondence among all the data sets with the help of the Iterative Closest Point method, and finally, we reconstruct the missing parts of the partial bones by creating a PCA model of the bone using the training data in correspondence.
</div></p>

More details on this project can be found in the attached report.
<img width="435" alt="Screenshot 2023-07-22 at 00 48 30" src="https://github.com/MarijaStojchevska/Probabilistic-Shape-Modelling/assets/18449614/9ec3fdb1-a277-49b5-a497-6ee5ac1e97c6">


<p><div align="center"><img width="567" src="https://github.com/MarijaStojchevska/Probabilistic-Shape-Modelling/assets/18449614/9ec3fdb1-a277-49b5-a497-6ee5ac1e97c6.png" (https://github.com/MarijaStojchevska/Probabilistic-Shape-Modelling/assets/18449614/9ec3fdb1-a277-49b5-a497-6ee5ac1e97c6.png)" > </div><div align="center"><i>Figure 1: Example recostruction of a Femur bone.</i></div></p>
