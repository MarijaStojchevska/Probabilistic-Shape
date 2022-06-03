package example

import java.awt.Color
import java.io.File

import scalismo.ui.api.ScalismoUI

import scalismo.io._
import scalismo.geometry.Landmark
import scalismo.geometry._3D;
import scalismo.geometry.EuclideanVector;
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.statisticalmodel.LowRankGaussianProcess;
import scalismo.mesh.TriangleMesh

import AlignmentAndCompletionTrials.TrainingAndorganising._
import AlignmentAndCompletionTrials.AlignmentTrials
import AlignmentAndCompletionTrials.AlignmentTrials.MeshLandmark
import AlignmentAndCompletionTrials.AlignmentTrials.MeshLandmarkArrays
import AlignmentAndCompletionTrials.ModelTrials;


object Tasks {

    def main( args: Array[ String ] ) {

   // ____________________________________Algorithm Control_______________________________________//
    		val useUI_intermediateResults 	: Boolean = true;		  // in decrising overpowering order
    				val useUI_forTraining				: Boolean = false;			//
						val useUI_PartialShapeAdmin	: Boolean = true;			//
    		val useUI_forResults 						: Boolean = true;			//
    		
		    // use iterative alignment ( General Procrustes Analysis ) by specifying a iteration limit > 0
    		val alignmentIterations 	: Int = 0;
    		val trainingIterations		: Int = 5;
    		val fittingIterations			: Int = 5; // during the projection phase

    		// if you want to use a pure Gaussian Process the whole time (from training to projection)
    		val use_raw_GP	: Boolean = false;

    		// specify the relative path where data is stored and shoudl be save to
    		val locTrainingDataMesh 	: String = "data/meshes";
				val locTrainingDataLM 		: String = "data/landmark";
				val locPartialShapesMesh	: String = "data/partials";
				val locPartialResults 		: String = "data/results/partial_";
		    val locIntermediaResults 	: String = "data/temps";
   // ____________________________________________________________________________________________//
   // ___________________________________________Notes____________________________________________//
	 //			training( mesh location , landmark location ) → model + reference                       //
	 //			getPartialData( partial mesh location , partial landmark location ) → partials          //
	 // 	  alignToPartialMesh( model , partials , reference , posteriori iterations )              //
	 // 	  store back results                                                                      //
	 // 	                   	                                                                      //
	 // 	  anything else is pure administrative ( control of UI , ... )                            //
   // ____________________________________________________________________________________________//

		    // expand temporary location with a result folder to store intermediate evaluation reults
		    val locIntermediaResultsEval = locIntermediaResults + "/results";
		    { // establishing presence of output & temp folders
		    	val intermediateFolders 		= new java.io.File( locIntermediaResults ).mkdir();
		    	val intermediateFoldersEval = new java.io.File( locIntermediaResultsEval ).mkdir();
		    	val results 								= new java.io.File( locPartialResults ).mkdir();
		    }
		    
		    // initialise Scalismo related library objects
		    scalismo.initialize()

		    // declare main variables
		    var trainedModel : PredictionModelFull = null;
		    var trainedGP : PredictionGaussianFull = null;
		    var reference = MeshLandmark( null , null );
		    var setOfIncompleteShapes = MeshLandmarkArrays( null , null );
		    var ui : ScalismoUI = null;
		    var resultsMeshes : Array[ TriangleMesh[ _3D ] ] = null
		    
		    if( useUI_intermediateResults ) { ui = ScalismoUI(); } // initialise UI here only if needed
		    	
		    /* Training *******************************************************************************/
		    /* ****multiple versions for training available *******************************************/
		    /* ****to switch method used in training, change call inside TrainingAndOrganising.scala***/
    		if( useUI_intermediateResults && useUI_forTraining ){ // training with UI ouputs

    				if( use_raw_GP ){ // if a raw GP is used through-out the whole calculation
		    				trainedGP = training_raw( locTrainingDataMesh , locTrainingDataLM , ui , 
																					alignmentIterations , trainingIterations , 
																					locIntermediaResults );
		    				reference = trainedGP.reference;
    				} else { // standard Mesh-Model is used
								trainedModel = training( locTrainingDataMesh , locTrainingDataLM , ui , 
																				 alignmentIterations , trainingIterations , 
																				 locIntermediaResults );
		    				reference = trainedModel.reference;
    				}
    				
    		} else { // only compute training

    				if( use_raw_GP ){ // if a raw GP is used through-out the whole calculation
		    				trainedGP = training_raw( locTrainingDataMesh , locTrainingDataLM, 
																					alignmentIterations , trainingIterations , 
																					locIntermediaResults );
		    				reference = trainedGP.reference;
    				} else { // standard Mesh-Model is used
								trainedModel = training( locTrainingDataMesh , locTrainingDataLM , 
																				 alignmentIterations , trainingIterations , 
																				 locIntermediaResults );
		    				reference = trainedModel.reference;
    				}
    				
    		}

		    /* Target Data*****************************************************************************/
				if( useUI_intermediateResults && useUI_PartialShapeAdmin ){ // load partials with UI output
				    setOfIncompleteShapes = getPartialData( locPartialShapesMesh , reference , ui , 
				    																				alignmentIterations , locIntermediaResults );
    		} else { // load partials silently
		    		setOfIncompleteShapes = getPartialData( locPartialShapesMesh , reference , 
		    																						alignmentIterations , locIntermediaResults );
    		}

		    /* Evaluate Targets by Model***************************************************************/
		    /* ****multiple versions for evaluation available******************************************/
		    /* ****to switch method used for evaluation, change call below*****************************/
		    if( use_raw_GP ){ // evaluation methods using raw GP instead of model
		    	
		    		resultsMeshes = ModelTrials.projectPartialMesh_RawGP_Iterative_CP( 
					    																	trainedGP.model , setOfIncompleteShapes , reference , 
					    																	fittingIterations , locIntermediaResultsEval );
		    	
		    } else { // evaluation mehtods using standard Statistical Mesh Model
		    	
	    			if( useUI_intermediateResults ){ // use UI for partial results and primed models
    					
	    					resultsMeshes = ModelTrials.projectPartialMesh_Iterative_CP_optim_UI(
										 												trainedModel.model , setOfIncompleteShapes , reference , 
										 												fittingIterations , ui , locIntermediaResultsEval );
  					
    					
    				} else { // silent evaluation
    					
	    					resultsMeshes = ModelTrials.projectPartialMesh_Iterative_CP_optim( 
	    																				trainedModel.model , setOfIncompleteShapes , 
	    																				reference , fittingIterations , 
	    																				locIntermediaResultsEval );
    					
    				}
	    			
		    }
				println("missing information found for partial shapes")

				if( !useUI_forResults ) return;

				/* display results ************************************************************************/
				if( useUI_forResults && !useUI_intermediateResults ){
						ui = ScalismoUI()
				}

				val setOfIncompleteShapesGroup = ui.createGroup("Partial Shapes")
				
				var counter : Int = 0;
				val setOfIncompleteShapesMeshesViews = setOfIncompleteShapes.meshes.map( mesh ⇒ {
				
						counter = counter + 1;
						ui.show( setOfIncompleteShapesGroup , mesh , "partial " + counter ) 
				
				} );    

				setOfIncompleteShapesMeshesViews.map( view ⇒ view.color =  new Color( 245 , 20 , 50 ) );

				counter = 0;
				val setOfIncompleteShapesLMsViews = setOfIncompleteShapes.lms.map( lm ⇒ {
						
						counter = counter + 1;
						ui.show( setOfIncompleteShapesGroup , lm , "partial " + counter + "L" ) 
													
				} );
				counter = 0;
				val resultsMeshesViews = resultsMeshes.map( mesh ⇒ {
				
						counter = counter + 1;
						ui.show( setOfIncompleteShapesGroup , mesh , "partial " + counter + "R" )

				} );
				resultsMeshesViews.map( view ⇒ view.color =  new Color( 100 , 200 , 70 ) );

				counter = 1;
				resultsMeshes.map( result ⇒ {
				
						MeshIO.writeMesh( result , 
															new java.io.File( locPartialResults + counter + ".stl" ) );
						counter = counter + 1;
					
				} );
				
    }

}