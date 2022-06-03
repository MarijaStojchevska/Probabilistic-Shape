package example

import java.awt.Color
import java.io.File

import breeze.linalg.{DenseMatrix, DenseVector}
import scalismo.common.{DiscreteField, Field, NearestNeighborInterpolator, PointId, RealSpace, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector, Landmark, Point, SquareMatrix, _3D}
import scalismo.io.{LandmarkIO, MeshIO}
import scalismo.kernels.{DiagonalKernel, GaussianKernel, MatrixValuedPDKernel, PDKernel}
import scalismo.mesh.{TriangleMesh, TriangleMesh3D}
import scalismo.numerics.PivotedCholesky.RelativeTolerance
import scalismo.numerics.{RandomMeshSampler3D, UniformMeshSampler3D}
import scalismo.registration.{LandmarkRegistration, Transformation, TranslationTransform}
import scalismo.statisticalmodel
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.statisticalmodel.{DiscreteLowRankGaussianProcess, GaussianProcess, MultivariateNormalDistribution, NDimensionalNormalDistribution, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI

import scala.io.StdIn.readLine


object ExampleApp_try {

  def main(args: Array[String]) {

    // setting a seed for the random generator to allow for reproducible results
    implicit val rng = scalismo.utils.Random( 42 )

    // required to initialize native libraries (VTK, HDF5 ..)
    scalismo.initialize()


    //Remove a group view by asking the order from the user
    def removeGroupView(groupView: scalismo.ui.api.Group): Unit = {
      groupView.remove()
    }

    // Your application code goes below here. Below is a dummy application that reads a mesh and displays it

    // create a visualization window
    val ui = ScalismoUI()

    // _____________________________MESHES_____________________________________//

    // read a mesh from file
    val referenceMesh = MeshIO.readMesh(new File("datasets/femur.stl" ) ).get
    //println(""+ referenceMesh.pointSet.points.toIndexedSeq.size)   //18879

    // display it
    val femur = ui.createGroup("Femur")
    val meshView = ui.show(femur, referenceMesh, "Reference" )
    // change its color

    val data = "Courses_SSM2016_Training_step2/step2/meshes/"
    /*
            for (i <- 1 until totalMeshes) yield {
             val instance = MeshIO.readMesh( meshFiles(i) ).get
             val interpretation = ui.show( originalData, instance, "mesh" + i )
              (instance, interpretation)
           }
    */

    val meshFiles = new File(data).listFiles()
    //meshFiles.foreach(file => println(s"${file}\n"))

    val numOfMeshes = meshFiles.size
    // println( "Total number of meshes: " + totalMeshes )
    // println(s"number of files is ${numOfMeshes}")

    val originalData = ui.createGroup("50_Bones")
    var i: Int = 1
    val (meshes, views) = meshFiles.map(meshFile => {
      val mesh = MeshIO.readMesh(meshFile).get
      //val meshView = ui.show(originalData, mesh, "mesh" + i)
      i = i + 1
      (mesh, meshView) // return a tuple of the mesh and the associated view
    }).unzip


    def center(sampleMesh : IndexedSeq[scalismo.geometry.Point[_3D]]): EuclideanVector[_3D]=
    {
      val vector = sampleMesh.map{point => point.toVector}
      val total_mesh_points = sampleMesh.size
      (vector.reduce((vec1, vec2) => vec1 + vec2)) / total_mesh_points
    }  //return the centroid (3D vector)


    //find centroid for the reference bone
    val pointsReference = referenceMesh.pointSet.points.toIndexedSeq
    val referenceCentroid = center(pointsReference)
    //change the position of the reference mesh according to the new centroid
    // val translation = TranslationTransform[_3D](referenceCentroid.*(-1))
    // val translatedReferenceMesh = referenceMesh.transform(translation)
    // ui.show(femur,translatedReferenceMesh, "Reference_prime")


    //find centroid for the given bones
    val pointsData = meshes.map(mesh => (mesh.pointSet.points.toIndexedSeq))
    val dataCentroids = pointsData.map(mesh => center(mesh))
    //making the centroid of the given bones equal to the origin
    val translations = dataCentroids.map(centVec => TranslationTransform[_3D](centVec.*(-1)))
    //change the position of the meshes according to their new centroids
    val translatedMeshes_center = (0 until numOfMeshes).map(i => meshes(i).transform(translations(i)))

    val points = translatedMeshes_center.map(mesh => (mesh.pointSet.points.toIndexedSeq))
    val centroids = points.map(mesh => center(mesh))
    //making the centroid of the given bones equal to the reference centroid
    val translation = centroids.map(cv => TranslationTransform[_3D](cv.+(referenceCentroid)))
    //change the position of the meshes according to their new centroids
    val translatedMeshes = (0 until numOfMeshes).map(i => translatedMeshes_center(i).transform(translation(i)))

    val movedData = ui.createGroup("50_Bones_translated")
    //(0 until numOfMeshes).map(i => ui.show(movedData, translatedMeshes(i), "mesh" + i))

    // _____________________________LANDMARKS_____________________________________//

    // val MeshPoints = meshes.map(mesh => (mesh.pointSet.points).toArray) //iterator to array

    //read and show the landmarks of the reference femur
    val landmarkReference = "datasets/femur.json"
    val landmarksFemur = LandmarkIO.readLandmarksJson[_3D](new java.io.File(landmarkReference)).get
    // val landmarkView = ui.show(femur, landmarksFemur, "landmarksFemur")

    //read and show the training data landmarks
    val landmarksAll = "Courses_SSM2016_Training_step2/step2/landmarks/"
    val landmarkFiles = new java.io.File(landmarksAll).listFiles
    //landmarkFiles.foreach(file => println(s"${file}\n"))


    // val originalLandmarks = ui.createGroup("50_Landmarks")
    var j: Int = 0
    val (landmarks) = landmarkFiles.map(landFile => {
      val landmark = LandmarkIO.readLandmarksJson[_3D](new java.io.File(landFile.getAbsolutePath)).get
      //val landmarkView = ui.show(originalData, landmark, "L" + j)
      j = j + 1
      (landmark) // return a tuple of the landmark and the associated view
    })


    //access the coordinates
    val coordinates =landmarks.map(land => land.toVector.map(coord => coord.point.toVector))
    //println(coordinates(0))
    val numCoordinates=coordinates.size  //50
    val numLMperMesh= coordinates(0).size //6 --> l0....l5


    //coordinates(landmarks_femur_1).(no_of_landmark).(no_of_coordinate_xyz)
    //using the same translations as for the meshes move the coordinates of the landmarks
    val LmPointsTranslation = (0 until numCoordinates).map(noMesh =>{
      (0 until numLMperMesh).map(noLM =>
        translation(noMesh)(translations(noMesh)(scalismo.geometry.Point(coordinates(noMesh)(noLM)(0),coordinates(noMesh)(noLM)(1),coordinates(noMesh)(noLM)(2)))))
    })

    //give IDs for each set of landmarks
    val LMs = LmPointsTranslation.map(meshLMs => {
      var n = -1
      meshLMs.map(point => {n=n+1; Landmark(s"L${n}", point)})
    })

    //show the translated landmarks
    // LMs.map(lms => {lms.map(point => ui.show(movedData,point,point.id))})

    // _____________________________ Rigid Alignment_____________________________________//

    //Rigid Transformation
    val rigidTransformation  = (0 until numOfMeshes).map (i =>
      LandmarkRegistration.rigid3DLandmarkRegistration(LMs(i), landmarksFemur, referenceCentroid.toPoint))

    //Landmarks alignment
    val alignedLM = (0 until numOfMeshes).map (i => {
      LMs(i).map(lm => lm.transform(rigidTransformation(i)))
    })

    //Show the rigidly transformed landmarks
    val alignedData = ui.createGroup("50_Bones_aligned")
    // val ref_LM_View = landmarksFemur.map(lm => ui.show(alignedData,lm,lm.id))
    //  val aligned_LM_Views = alignedLM.map(eachSet => {eachSet.map(lm => ui.show(alignedData,lm,lm.id))})

    //Meshes alignment
    val alignedMeshes = (0 until numOfMeshes).map(i =>
      translatedMeshes(i).transform(rigidTransformation(i)))  //here we use the translated meshes

    //for faster and easier rigid transformation
    //Show the rigidly transformed meshes
    // val ref_Mesh_View = ui.show(alignedData, referenceMesh, "reference_mesh")
    //  val aligned_Meshes_Views = (0 until numOfMeshes).map(i => ui.show(alignedData, alignedMeshes(i), "mesh" + i))

    // _____________________________Save the aligned bones and landmarks in files_____________________________________//
    /*  val meshPath = "./Aligned_dataset/meshes/"
      val lmPath = "./Aligned_dataset/landmarks/"

      MeshIO.writeMesh(referenceMesh, new java.io.File(meshPath + "referenceMesh.stl"))
      (0 until numOfMeshes).foreach(i => MeshIO.writeMesh((alignedMeshes(i)), new File(meshPath + s"${i}.stl")))

      LandmarkIO.writeLandmarksJson[_3D](landmarksFemur, new java.io.File(lmPath + "referenceLM.json"))
      (0 until numOfMeshes).foreach(i => LandmarkIO.writeLandmarksJson[_3D](alignedLM(i), new File(lmPath + s"${i}.json")))
  */



    // __________________1) Create correspondence between the training data through using ICP__________________________//

    removeGroupView(femur)
    //Define mean and covariance function (kernel) for the Gaussian Process

    //Mean
    val zeroMean = scalismo.common.VectorField(RealSpace[_3D], (pt : scalismo.geometry.Point[_3D]) => EuclideanVector.zeros[_3D])

    val scalarValuedKernel = GaussianKernel[_3D](25) * 50

    case class XmirroredKernel(ker : PDKernel[_3D]) extends PDKernel[_3D] {
      override def domain : RealSpace[_3D] = RealSpace[_3D]
      override def k(x: Point[_3D], y: Point[_3D]) = ker(Point(x(0) * -1f ,x(1), x(2)), y)
    }
    def SymmetrizeKernel(ker : PDKernel[_3D]) : MatrixValuedPDKernel[_3D] = {
      val xmirrored = XmirroredKernel(ker)

      val k2 = DiagonalKernel(xmirrored * -1f, xmirrored, xmirrored)
      k2
    }
    val sim = SymmetrizeKernel(scalarValuedKernel)





   // val diagKernel1 = DiagonalKernel[_3D](GaussianKernel(sigma = 30)*10, outputDim = 3)  //type two

    //type three
    /*
        case class ChangePointKernel(kernel1 : MatrixValuedPDKernel[_3D], kernel2 : MatrixValuedPDKernel[_3D])
          extends MatrixValuedPDKernel[_3D]() {

          override def domain = RealSpace[_3D]
          val outputDim = 3

          def s(p: Point[_3D]) =  1.0 / (1.0 + math.exp(-p(0)))
          def k(x: Point[_3D], y: Point[_3D]) = {
            val sx = s(x)
            val sy = s(y)
            kernel1(x,y) * sx * sy + kernel2(x,y) * (1-sx) * (1-sy)
          }

        }

        val gk1 = DiagonalKernel(GaussianKernel[_3D](140.0), 3)
        val gk2 = DiagonalKernel(GaussianKernel[_3D](150.0), 3)
        val gk3 = DiagonalKernel(GaussianKernel[_3D](160.0), 3)
        val gk4 = DiagonalKernel(GaussianKernel[_3D](170.0), 3)
        val a = ChangePointKernel(gk3, gk4)
        val b = ChangePointKernel(gk1, gk2)
        val sum = a + b

     */

    val gp = GaussianProcess(zeroMean, sim)


    val firstAligned = alignedMeshes(0)


    val sampler = RandomMeshSampler3D(
      firstAligned,
      numberOfPoints = 300,
      seed = 45)


    val lowRankGP = scalismo.statisticalmodel.LowRankGaussianProcess.approximateGPNystrom(gp, sampler, numBasisFunctions = 100) //finite rank decomposition of the GP, retaining only the 100 most prominent basis functions

    val deformedGroup = ui.createGroup("deformed")
    val  model = StatisticalMeshModel(firstAligned, lowRankGP) //OUR MODEL
    val reference = model.mean //OUR REFERENCE


    val lm_ref_points = reference.pointSet
    val lm_ref_ids = lm_ref_points.pointIds.toIndexedSeq

    /*
    val referencePoints = UniformMeshSampler3D(reference, 1000).sample.map(s => s._1)
    val lm_ref_ids = referencePoints.map(lm => reference.pointSet.findClosestPoint(lm).id) //identifiers for the blue points
    val lm_ref_points = lm_ref_ids.map(id => reference.pointSet.point(id))
     */

    ui.show(deformedGroup, model, "MODEL")




    def attributeCorrespondences(currentMesh : TriangleMesh3D, target : TriangleMesh3D) : Seq[(PointId, scalismo.geometry.Point[_3D])]  = {

      val correspond = lm_ref_ids.map{id =>
        val pt = currentMesh.pointSet.point(id)
        val clostestPtsOnTarget = target.pointSet.findClosestPoint(pt).point
        (id, clostestPtsOnTarget)
      }

      val newPoints =  lm_ref_ids.map(id => currentMesh.pointSet.point(id))
      //ui.show(MovingMeshPoints, newPoints, "MovingMeshPoints")

      correspond
    }

    //finding and showing the candidate correspondences for each of the target meshes
    val correspondences = ui.createGroup("Correspondences")

    def CandidateCorrespondences(currentMesh : TriangleMesh3D, targetNo : Int) :  Seq[(PointId, scalismo.geometry.Point[_3D])]  = {
      val candidates = attributeCorrespondences(currentMesh, alignedMeshes(targetNo))

      val newCandidates =  candidates.map(c => c._2).toIndexedSeq
      //val candidatePoints = ui.show(correspondences, newCandidates, "newCorrespondences" + targetNo)
      //candidatePoints.color =  new Color(245,20,50)
      candidates
    }

    //create deformation field between landmarks of the reference <----> candidate points of the target mesh (for each target mesh)

    //CHANGED
    def deformation(getCandidates : scala.Seq[scalismo.geometry.Point[_3D]], dataPoints : scala.Seq[scalismo.geometry.Point[_3D]]) : IndexedSeq[EuclideanVector[_3D]] ={
      val candidates = getCandidates
      val deform = (0 until lm_ref_points.points.size).map(j => (candidates(j) - dataPoints(j)))
      deform
    }

    val deformationView = ui.createGroup("Deformation Fields")


    //CHANGED
    def deformationFields(getCandidates : scala.Seq[scalismo.geometry.Point[_3D]], currentMesh : TriangleMesh3D) : DiscreteField[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]]={

      val newPoints = lm_ref_ids.map{id => currentMesh.pointSet.point(id)}
      val values =  deformation(getCandidates, newPoints)

      val deformFields = DiscreteField[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]](UnstructuredPointsDomain(newPoints), values.toIndexedSeq)

      ui.show(deformationView, deformFields, "deformation")

      deformFields
    }

    //___________________________Using the candidate deformation fields for regression________________________________//

    def trainingData(candidateCorresp : Seq[(PointId, scalismo.geometry.Point[_3D])],  noiseParameter : Double):  Seq[(PointId, scalismo.geometry.Point[_3D], MultivariateNormalDistribution)]= {

      val noise = MultivariateNormalDistribution(DenseVector.zeros[Double](3),  noiseParameter * DenseMatrix.eye[Double](3))
      val dataGenerator = candidateCorresp.map(c => (c._1, c._2, noise))
      dataGenerator
    }


    def posteriorFunction(candidateCorresp : Seq[(PointId, scalismo.geometry.Point[_3D])],  noiseParameter : Double, currentMesh : TriangleMesh3D, noIter : Int, iter : Int) : TriangleMesh[_3D] = {

      val trainData = trainingData(candidateCorresp, noiseParameter) //CHANGED, use  lm_ref_points

      val posteriorGP = model.posterior(trainData.toIndexedSeq)

      val newMean = posteriorGP.mean


      val posteriorGroup = ui.createGroup("Posterior")
      ui.show(posteriorGroup, newMean, "Mean Deformation")
      if(iter < noIter)
      {
        // Thread.sleep(3000)
        removeGroupView(posteriorGroup)
      }
      else
      {
        ui.show(posteriorGroup, posteriorGP, "StatisticalModel")
      }


      //  ui.show(posteriorGroup, currentMesh, "Mesh Defrmation" )
      //  ui.addTransformation(posteriorGroup, newMean, "transform")


      newMean
    }




    val dataPath= "./Data_ICP/data/"
    //___________________________Recursion for Iteration________________________________//
    def recursion(currentMesh : TriangleMesh[_3D], nbIterations : Int, targetNo : Int  , iter : Int, noiseParameter : Double):Any= {

      if(targetNo >= alignedMeshes.size)  // change partials.size to 1 for TESTING
      {
        println("END RECURSION")
        return
      }

      //compute candidates
      val candidates = CandidateCorrespondences(currentMesh, targetNo) //returnes ids +  new candidates on the target Mesh!

      //compute deformation fields
      val newCandidates =  candidates.map(c => c._2).toIndexedSeq
      // deformationFields(newCandidates, currentMesh)

      //compute posterior mean
      val mean = posteriorFunction(candidates, noiseParameter, currentMesh, nbIterations, iter)



      val avgDistance  = scalismo.mesh.MeshMetrics.avgDistance(alignedMeshes(targetNo), mean)
      val hausdorffDistance  = scalismo.mesh.MeshMetrics.hausdorffDistance(alignedMeshes(targetNo), mean)

      println("Number of iteration: "+ iter + ", for target mesh: " + targetNo)
      println("Average distance: "+ avgDistance)
      println("Hausdorff distance: " + hausdorffDistance)
      println(" ")

      if(iter < nbIterations)
      {
        recursion(mean, nbIterations, targetNo, iter + 1, noiseParameter) //reducing the noise parameter gives us slightly better results
      }
      else{
        MeshIO.writeMesh(mean, new java.io.File(dataPath + s"mean${targetNo}.stl"))
        recursion(reference, nbIterations, targetNo+1, 0, noiseParameter) //CHANGED
      }

    }

    println("START RECURSION")

    recursion(reference, 1, 0, 0,  1.0) ///CHANGED




        // __________________2) Use the training data in correspondence to create a training model by using PCA__________________________//

            val pathNewTraining = "./Data_ICP/data2"
            val readNewTraining = new File(pathNewTraining).listFiles()


            val generatedData = readNewTraining.map(meshFile => {MeshIO.readMesh(meshFile).get})

            val firstInstance = generatedData(0)
            val firstInstancePoints =  firstInstance.pointSet

            val deFields = generatedData.map{each =>
              val vectors = firstInstance.pointSet.pointIds.map{id : PointId   =>
                each.pointSet.point(id) - firstInstance.pointSet.point(id)
              }.toIndexedSeq
              DiscreteField[_3D, UnstructuredPointsDomain[_3D], EuclideanVector[_3D]](firstInstancePoints, vectors)
            }

            val interpolator = NearestNeighborInterpolator[_3D, EuclideanVector[_3D]]()
            val continuousFields  = deFields.map(f => f.interpolate(interpolator))

            val GP = DiscreteLowRankGaussianProcess.createUsingPCA(firstInstancePoints, continuousFields , RelativeTolerance(1e-8))


            val pcaModel = StatisticalMeshModel(firstInstance, GP.interpolate(interpolator))
            val pcaGroup = ui.createGroup("PCA")
            ui.show(pcaGroup, pcaModel, "PCA Model")



        // __________________3) Use the training model to reconstruct the missing parts__________________________//


        // _____________________________Partials Reading_____________________________________//

        //reading and showing the partial meshes
        val partialsPath = "Partials/meshes"
        val partialsFiles = new File(partialsPath).listFiles()

        val partialMeshes = ui.createGroup("10_Partials")
        var p: Int = 0
        val partials = partialsFiles.map(partFile => {
          val partial = MeshIO.readMesh(partFile).get
          ui.show(partialMeshes, partial, "partial" + p)
          p = p + 1
          partial
        })

        println("partials Done")




































  }
}
