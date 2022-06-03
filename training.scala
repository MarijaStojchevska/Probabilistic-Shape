package example

import java.awt.Color

import scalismo.ui.api.ScalismoUI

import scalismo.io._
import scalismo.geometry.Landmark
import scalismo.geometry._3D;
import scalismo.geometry.EuclideanVector;
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.statisticalmodel.LowRankGaussianProcess;

import AlignmentAndCompletionTrials.AlignmentTrials
import AlignmentAndCompletionTrials.AlignmentTrials.MeshLandmark
import AlignmentAndCompletionTrials.AlignmentTrials.MeshLandmarkArrays
import AlignmentAndCompletionTrials.ModelTrials;
import scalismo.mesh.TriangleMesh


package AlignmentAndCompletionTrials{

object TrainingAndorganising {

		case class PredictionModelFull( model 		: StatisticalMeshModel , 
														 				reference	: MeshLandmark ){}
		
		case class PredictionGaussianFull( model 		  : LowRankGaussianProcess[ _3D , EuclideanVector[ _3D ] ] , 
														 					 reference	: MeshLandmark ){}
	
		
		/**
		 * by extracting the training elements out of the remaining calculation sequence, we can
		 * free active memory space for our calculations
		 */
		def training( meshLocations : String , landmarkLocations : String , alignmentIterations : Int , 
									modellingIterations : Int , tempFolder : String ) : PredictionModelFull = {
				
				// declaring main variables
				var trainingSet 		: MeshLandmarkArrays 	= null;
				var reference 			: MeshLandmark 					= null;
				var predictionModel : StatisticalMeshModel 	= null;

				// expanding temporary locations for intermediate steps
		    val locIntermediaAlignedTrainingM 	: String = tempFolder + "/alignedTrain/meshes";
		    val locIntermediaAlignedTrainingL 	: String = tempFolder + "/alignedTrain/LMs";
				val locIntermediaModel						 	: String = tempFolder + "/model.h5";
				val locIntermediaRefLM						 	: String = tempFolder + "/reference.json";
				val locIntermediaRefM						 	  : String = tempFolder + "/reference.stl";
		    { // establishing presence of temp folders
		    	val intermediateFoldersM = new java.io.File( locIntermediaAlignedTrainingM ).mkdir();
		    	val intermediateFoldersL = new java.io.File( locIntermediaAlignedTrainingL ).mkdir();
		    }
		    
		    { // read precomputed model
				    val model = ModelTrials.readModelFile( locIntermediaModel );
		    		val ref		= AlignmentTrials.readMeshLM( locIntermediaRefM , locIntermediaRefLM );
		    		// if pre-computed model could be read, we are done
		    		if( model != null && ref != null ) return PredictionModelFull( model , ref );
		    }

		    // read precomputed aligned data
				trainingSet = AlignmentTrials.readMeshsLMs( locIntermediaAlignedTrainingM , 
																										locIntermediaAlignedTrainingL );
					
				/* reading and preparing training set (aligning meshes and defining reference,...)*********/
				if( trainingSet.meshes == null ){ // no precomputed aligned data found
					
						trainingSet = AlignmentTrials.readMeshsLMs( "data/meshes" , "data/landmarks" );
				    println("reading of training data finished")
		
				    /*
				     *  this centroid alignment uses all shape points instead of only landmarks
				     *  this step speeds up a the following Procrustes Analysis (alignment)
				     */
				    AlignmentTrials.centeroidAlignment( trainingSet );
				    println("centeroid alignment of training data finished")
				    
				    if( alignmentIterations > 0 ){ // valid number of iterations : choose iterative alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks( trainingSet , 
				    																												 alignmentIterations );
				    } else { // non-valid number of iterations : choose single-iteration alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks_singleRun( trainingSet );
				    }
				    println("aligning of training data finished")
				    
				    AlignmentTrials.writeMeshsLMs( locIntermediaAlignedTrainingM , 
				    															 locIntermediaAlignedTrainingL , trainingSet );
					
				}
		    
				/* training a model based on the training data*********************************************/
		    /* ****multiple versions available*********************************************************/
		    /* ****to switch method, change call below*************************************************/
		    predictionModel = ModelTrials. buildGaussProcessModel_ICP_relaxed_PCA2( 
		    																					trainingSet , reference , modellingIterations );

		    println("prediction model created")
		   	ModelTrials.writeModelFile( predictionModel , locIntermediaModel );

				return PredictionModelFull( predictionModel , reference );
			
		}
		
		/**
		 * by extracting the training elements out of the remaining calculation sequence, we can
		 * free active memory space for our calculations
		 * the raw version returns a pure GP
		 */
		def training_raw( meshLocations : String , landmarkLocations : String , 
											alignmentIterations : Int , modellingIterations : Int , tempFolder : String ) : 
																																					PredictionGaussianFull = {

				// declaring main variables
				var trainingSet 	: MeshLandmarkArrays 																		 = null;
				var reference 		: MeshLandmark 																					 = null;
				var predictionGP 	: LowRankGaussianProcess[ _3D , EuclideanVector[ _3D ] ] = null;
			
				// expanding temporary locations for intermediate steps
		    val locIntermediaAlignedTrainingM 	: String = tempFolder + "/alignedTrain/meshes";
		    val locIntermediaAlignedTrainingL 	: String = tempFolder + "/alignedTrain/LMs";
				val locIntermediaModel						 	: String = tempFolder + "/model.h5";
				val locIntermediaRefLM						 	: String = tempFolder + "/reference.json";
				val locIntermediaRefM						 	  : String = tempFolder + "/reference.stl";
		    { // establishing presence of temp folders
		    	val intermediateFoldersM = new java.io.File( locIntermediaAlignedTrainingM ).mkdir();
		    	val intermediateFoldersL = new java.io.File( locIntermediaAlignedTrainingL ).mkdir();
		    }
		    
		    { // read precomputed model
		    		// error : problem reading raw Gaussian Process from file, we get different repres.
		    		val gp : LowRankGaussianProcess[ _3D , EuclideanVector[ _3D ] ] = null;
//				    																		ModelTrials.readModelFile( locIntermediaModel ).gp;
		    		val ref		= AlignmentTrials.readMeshLM( locIntermediaRefM , locIntermediaRefLM );
		    		// if pre-computed model could be read
		    		if( gp != null && ref != null ) return PredictionGaussianFull( gp , ref );
		    }

		    // read precomputed aligned data
				trainingSet = AlignmentTrials.readMeshsLMs( locIntermediaAlignedTrainingM , 
																										locIntermediaAlignedTrainingL );
																											
				/* reading and preparing training set (aligning meshes and defining reference,...)*********/
				if( trainingSet.meshes == null ){ // no precomputed aligned data found
			
						var trainingSet = AlignmentTrials.readMeshsLMs( "data/meshes" , "data/landmarks" );
				    println("reading of training data finished")
		
				    /*
				     *  this centroid alignment uses all shape points instead of only landmarks
				     *  this step speeds up a the following Procrustes Analysis (alignment)
				     */
				    AlignmentTrials.centeroidAlignment( trainingSet );
				    println("centeroid alignment of training data finished")
				    
						var reference : MeshLandmark = null;
				    if( alignmentIterations > 0 ){ // valid number of iterations : choose iterative alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks( trainingSet , 
				    																												 alignmentIterations );
				    } else { // non-valid number of iterations : choose single-iteration alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks_singleRun( trainingSet );		    	
				    }
				    println("aligning of training data finished")
	
				    AlignmentTrials.writeMeshsLMs( locIntermediaAlignedTrainingM , 
				    															 locIntermediaAlignedTrainingL , trainingSet );

				}

				/* training a GP based on the training data************************************************/
		    /* ****multiple versions available*********************************************************/
		    /* ****to switch method, change call below*************************************************/
				predictionGP = ModelTrials.createGaussianProcess_HighFlexKernel( reference , 30 , 35 , false );
		    println("prediction model created")
				
		    // to use predefined IO methods, GP has to be expanded to Statistical Mesh Model
		    //  --> storing data in associated wrapper
		    val predictionModel = PredictionModelFull( 
		    													StatisticalMeshModel( reference.mesh , predictionGP ) , reference )
		   	ModelTrials.writeModelFile( predictionModel.model , locIntermediaModel );

				return PredictionGaussianFull( predictionGP , reference );
			
		};

		/**
		 * by extracting the training elements out of the remaining calculation sequence, we can
		 * free active memory space for our calculations
		 * this UI version displays intermediate results
		 */
		def training( meshLocations : String , landmarkLocations : String , UI : ⇒ ScalismoUI , 
									alignmentIterations : Int , modellingIterations : Int , tempFolder : String ) : 
																																						PredictionModelFull = {
													
				val trainingGroup = UI.createGroup("training Shapes and Results")

				// declaring main variables
				var trainingSet 		: MeshLandmarkArrays 		= null;
				var reference 			: MeshLandmark 					= null;
				var predictionModel : StatisticalMeshModel 	= null;
			
				// expanding temporary locations for intermediate steps
		    val locIntermediaAlignedTrainingM 	: String = tempFolder + "/alignedTrain/meshes";
		    val locIntermediaAlignedTrainingL 	: String = tempFolder + "/alignedTrain/LMs";
				val locIntermediaModel						 	: String = tempFolder + "/model.h5";
				val locIntermediaRefLM						 	: String = tempFolder + "/reference.json";
				val locIntermediaRefM						 	  : String = tempFolder + "/reference.stl";
		    { // establishing presence of temp folders
		    	val intermediateFoldersM = new java.io.File( locIntermediaAlignedTrainingM ).mkdir();
		    	val intermediateFoldersL = new java.io.File( locIntermediaAlignedTrainingL ).mkdir();
		    }
		    
		    { // read precomputed model
				    val model = ModelTrials.readModelFile( locIntermediaModel );
		    		val ref		= AlignmentTrials.readMeshLM( locIntermediaRefM , locIntermediaRefLM );
		    		if( model != null && ref != null ) { // if pre-computed model could be read
		    			
								UI.show( trainingGroup , predictionModel , "sampler") 
		    				return PredictionModelFull( model , ref );
		    		
		    		}
		    }

		    // read precomputed aligned data
				trainingSet = AlignmentTrials.readMeshsLMs( locIntermediaAlignedTrainingM , 
																										locIntermediaAlignedTrainingL );
																											
				/* reading and preparing training set (aligning meshes and defining reference,...)*********/
				if( trainingSet.meshes == null ){ // no precomputed aligned data found
				
						trainingSet = AlignmentTrials.readMeshsLMs( "data/meshes" , "data/landmarks" );
				    println("reading of training data finished")
		
				    var counter : Int = 0;
						val trainingMeshesViews1 = trainingSet.meshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( trainingGroup , mesh , "partial " + counter + " original" ) 
						
						} );
		
				    /*
				     *  this centroid alignment uses all shape points instead of only landmarks
				     *  this step speeds up a the following Procrustes Analysis (alignment)
				     */
				    AlignmentTrials.centeroidAlignment( trainingSet );
				    println("centeroid alignment of training data finished")
		
				    counter = 0;
						val trainingMeshesViews2 = trainingSet.meshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( trainingGroup , mesh , "partial " + counter + " CentAlign" ) 
						
						} );
		
				    if( alignmentIterations > 0 ){ // valid number of iterations : choose iterative alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks( trainingSet , alignmentIterations );
				    } else { // non-valid number of iterations : choose single-iteration alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks_singleRun( trainingSet );		    	
				    }
				    println("aligning of training data finished")

				    AlignmentTrials.writeMeshsLMs( locIntermediaAlignedTrainingM , 
				    															 locIntermediaAlignedTrainingL , trainingSet );

				    counter = 0;
						val trainingMeshesViews3 = trainingSet.meshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( trainingGroup , mesh , "partial " + counter + " rigid" ) 
						
						} );

				}

				/* training a model based on the training data*********************************************/
		    /* ****multiple versions available*********************************************************/
		    /* ****to switch method, change call below*************************************************/
				predictionModel = ModelTrials.buildGaussProcessModel_ICP_relaxed_PCA2( 
																									trainingSet , reference , modellingIterations )
		    println("prediction model created")

		   	ModelTrials.writeModelFile( predictionModel , locIntermediaModel );

				UI.show( trainingGroup , predictionModel , "sampler") 

				return PredictionModelFull( predictionModel , reference );
			
		};

		/**
		 * by extracting the training elements out of the remaining calculation sequence, we can
		 * free active memory space for our calculations
		 * the raw version returns a pure GP
		 * this UI version displays intermediate results
		 */
		def training_raw( meshLocations : String , landmarkLocations : String , UI : ⇒ ScalismoUI , 
											alignmentIterations : Int , modellingIterations : Int , tempFolder : String ) : 
																																					PredictionGaussianFull = {
													
				val trainingGroup = UI.createGroup("training Shapes and Results")

				// declaring main variables
				var trainingSet : MeshLandmarkArrays = null;
				var reference : MeshLandmark = null;
				var predictionGP : LowRankGaussianProcess[ _3D , EuclideanVector[ _3D ] ] = null;
			
				// expanding temporary locations for intermediate steps
		    val locIntermediaAlignedTrainingM 	: String = tempFolder + "/alignedTrain/meshes";
		    val locIntermediaAlignedTrainingL 	: String = tempFolder + "/alignedTrain/LMs";
				val locIntermediaModel						 	: String = tempFolder + "/model.h5";
				val locIntermediaRefLM						 	: String = tempFolder + "/reference.json";
				val locIntermediaRefM						 	  : String = tempFolder + "/reference.stl";
		    { // establishing presence of temp folders
		    	val intermediateFoldersM = new java.io.File( locIntermediaAlignedTrainingM ).mkdir();
		    	val intermediateFoldersL = new java.io.File( locIntermediaAlignedTrainingL ).mkdir();
		    }
		    
		    { // read precomputed model
		    		// error : problem reading raw Gaussian Process from file, we get different repres.
				    val gp : LowRankGaussianProcess[ _3D , EuclideanVector[ _3D ] ] = null;
//				    																			ModelTrials.readModelFile( locIntermediaModel ).gp;
		    		val ref		= AlignmentTrials.readMeshLM( locIntermediaRefM , locIntermediaRefLM );
		    		if( gp != null && ref != null ) { // if pre-computed model could be read
		    			
							UI.show( trainingGroup , predictionGP , "sampler") 
		    			return PredictionGaussianFull( gp , ref );
		    		
		    		}
		    }

		    // read precomputed aligned data
				trainingSet = AlignmentTrials.readMeshsLMs( locIntermediaAlignedTrainingM , 
																										locIntermediaAlignedTrainingL );
																											
				/* reading and preparing training set (aligning meshes and defining reference,...)*********/
				if( trainingSet.meshes == null ){ // no precomputed aligned data found
			
						trainingSet = AlignmentTrials.readMeshsLMs( "data/meshes" , "data/landmarks" );
				    println("reading of training data finished")
		
				    var counter : Int = 0;
						val trainingMeshesViews1 = trainingSet.meshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( trainingGroup , mesh , "partial " + counter + " original" ) 
						
						} );
		
				    /*
				     *  this centroid alignment uses all shape points instead of only landmarks
				     *  this step speeds up a the following Procrustes Analysis (alignment)
				     */
				    AlignmentTrials.centeroidAlignment( trainingSet );
				    println("centeroid alignment of training data finished")
		
				    counter = 0;
						val trainingMeshesViews2 = trainingSet.meshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( trainingGroup , mesh , "partial " + counter + " CentAlign" ) 
						
						} );
		
				    if( alignmentIterations > 0 ){// valid number of iterations : choose iterative alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks( trainingSet , alignmentIterations );
				    } else {// non-valid number of iterations : choose single-iteration alignment
				    		reference = AlignmentTrials.alignMeshesAndLandmarks_singleRun( trainingSet );		    	
				    }
				    println("aligning of training data finished")

				    AlignmentTrials.writeMeshsLMs( locIntermediaAlignedTrainingM , 
															 locIntermediaAlignedTrainingL , trainingSet );

				    counter = 0;
						val trainingMeshesViews3 = trainingSet.meshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( trainingGroup , mesh , "partial " + counter + " rigid" ) 
						
						} );
						
				}

				/* training a model based on the training data*********************************************/
		    /* ****multiple versions available*********************************************************/
		    /* ****to switch method, change call below*************************************************/
				predictionGP = ModelTrials.createGaussianProcess_HighFlexKernel( reference , 30 , 35 , false );
		    println("prediction model created")

		    val predictionModel = PredictionModelFull( 
		    													StatisticalMeshModel( reference.mesh , predictionGP ) , reference )
		   	ModelTrials.writeModelFile( predictionModel.model , locIntermediaModel );

				UI.show( trainingGroup , predictionGP , "sampler") 

				return PredictionGaussianFull( predictionGP , reference );
			
		}

		/**
		 * by extracting the training elements out of the remaining calculation sequence, we can
		 * free active memory space for our calculations
		 */
		def getPartialData( location : String , reference : MeshLandmark , alignmentIterations : Int , 
												tempLocation : String ) : MeshLandmarkArrays = {

				// declaring main variables
				var setOfIncompleteShapes : MeshLandmarkArrays = null;								
			
				// expanding temporary locations for intermediate steps
				val locIntermediaAlignedPartialsM 	: String = tempLocation + "/alignedPartial/meshes";
				val locIntermediaAlignedPartialsL 	: String = tempLocation + "/alignedPartial/landmarks";
		    { // establishing presence of temp folders
		    	val intermediateFoldersM = new java.io.File( locIntermediaAlignedPartialsM ).mkdir();
		    	val intermediateFoldersL = new java.io.File( locIntermediaAlignedPartialsL ).mkdir();
		    }

		    // read precomputed aligned data
				setOfIncompleteShapes = AlignmentTrials.readMeshsLMs( locIntermediaAlignedPartialsM , 
																															locIntermediaAlignedPartialsL );
																					
				// if precomputed aligned data found
				if( setOfIncompleteShapes.meshes != null ) return setOfIncompleteShapes;
			
				/* reading & prepare partial set (aligning meshes and defining reference,...)**************/
				var setOfIncompleteShapesMeshes = AlignmentTrials.readMeshes( location );
		    println("reading of partial shape data finished")
		    
		    // generate a suiteable set of landmarks based on the reference for the partial data
				var setOfIncompleteShapesLMs    = AlignmentTrials.alignLandmarksTo( reference.lm , 
																																			setOfIncompleteShapesMeshes );
		    
		    // store partial information in wrapper
		    setOfIncompleteShapes = MeshLandmarkArrays( setOfIncompleteShapesMeshes , 
		    																						setOfIncompleteShapesLMs )

				// align partials to reference / trained model
		  	AlignmentTrials.alignMeshesAndLandmarksToReference( setOfIncompleteShapes , reference )
		    																								
		    println("aligning of partial shape data finished")
		    
		    AlignmentTrials.writeMeshsLMs( locIntermediaAlignedPartialsM , locIntermediaAlignedPartialsL , 
		    															 setOfIncompleteShapes );
		    return setOfIncompleteShapes;
			
		}

		/**
		 * by extracting the training elements out of the remaining calculation sequence, we can
		 * free active memory space for our calculations
		 * this UI version displays intermediate results
		 */
		def getPartialData( location : String , reference : MeshLandmark , UI : ⇒ ScalismoUI , 
												alignmentIterations : Int , tempLocation : String ) : MeshLandmarkArrays = {
			
				val partialGroup = UI.createGroup("target Data")
				
				// declaring main variables
				var setOfIncompleteShapes : MeshLandmarkArrays = null;
				var counter								: Int								 = 0;
			
				// expanding temporary locations for intermediate steps
				val locIntermediaAlignedPartialsM 	: String = tempLocation + "/alignedPartial/meshes";
				val locIntermediaAlignedPartialsL 	: String = tempLocation + "/alignedPartial/landmarks";
		    { // establishing presence of temp folders
		    	val intermediateFoldersM = new java.io.File( locIntermediaAlignedPartialsM ).mkdir();
		    	val intermediateFoldersL = new java.io.File( locIntermediaAlignedPartialsL ).mkdir();
		    }

		    // read precomputed aligned data
				setOfIncompleteShapes = AlignmentTrials.readMeshsLMs( locIntermediaAlignedPartialsM , 
																															locIntermediaAlignedPartialsL );
																					
				/* reading & prepare partial set (aligning meshes and defining reference,...)**************/
				if( setOfIncompleteShapes.meshes != null ) {// no precomputed aligned data found
			
						var setOfIncompleteShapesMeshes = AlignmentTrials.readMeshes( location );
				    println("reading of partial shape data finished")
				    	
				    counter = 0;
						val partialMeshesViews1 = setOfIncompleteShapesMeshes.map( mesh ⇒ {
						
								counter = counter + 1;
								UI.show( partialGroup , mesh , "partial " + counter + " original" ) 
						
						} );
		
				    // generate a suiteable set of landmarks based on the reference for the partial data
						var setOfIncompleteShapesLMs    = AlignmentTrials.alignLandmarksTo( reference.lm , 
																																					setOfIncompleteShapesMeshes );
				    
				    // store partial information in wrapper
				    var setofIncompleteShapes = MeshLandmarkArrays( setOfIncompleteShapesMeshes , 
				    																								setOfIncompleteShapesLMs )
		
						// align partials to reference / trained model
				  	AlignmentTrials.alignMeshesAndLandmarksToReference( setofIncompleteShapes , reference )
		
				    AlignmentTrials.writeMeshsLMs( locIntermediaAlignedPartialsM , locIntermediaAlignedPartialsL , 
				    															 setOfIncompleteShapes );
				    
				}
		  	counter = 0;
				val partialMeshesViews2 = setOfIncompleteShapes.meshes.map( mesh ⇒ {
				
						counter = counter + 1;
						UI.show( partialGroup , mesh , "partial " + counter + " rigid" ) 
				
				} );

		    println("aligning of partial shape data finished")
		    
		    return setOfIncompleteShapes;
			
		}

}

}