package example;

import java.io._
import java.io.File

import scalismo.ui.api.ScalismoUI
import scalismo.ui.api.{ Group ⇒ apiGroup }

import scalismo.io.StatisticalModelIO;
import scalismo.registration._;
import scalismo.geometry._;
import Array._;
import scalismo.mesh._;
import scalismo.common._;
import scalismo.statisticalmodel._
import scalismo.registration._
import scalismo.statisticalmodel._
import scalismo.statisticalmodel.dataset._
import scalismo.numerics.UniformMeshSampler3D;
import scalismo.numerics.PivotedCholesky.RelativeTolerance
import scalismo.numerics.RandomMeshSampler3D
import scalismo.kernels._;
import breeze.linalg.{ DenseMatrix , DenseVector };
//import scala.collection.parallel.ParIterableLike.Reduce

import AlignmentAndCompletionTrials.AlignmentTrials._;

package AlignmentAndCompletionTrials {
	
object ModelTrials {
		
		/**
		 * read Statistical Mesh Model from file with the given file-path
		 * ( .h5 files native to Scalismo are presumed )
		 * returns null, if file does not exist !
		 */
		def readModelFile( name : String ) : StatisticalMeshModel = {
			
			val file = new java.io.File( name );
			if( !file.exists() ) return null;

			val model = StatisticalModelIO.readStatisticalMeshModel( file ).get
			
			return model;
			
		}

		/**
		 * write Statistical Mesh Model to file with the given file-path
		 * ( .h5 files native to Scalismo are presumed )
		 */
		def writeModelFile( model : StatisticalMeshModel , name : String ) {

			StatisticalModelIO.writeStatisticalMeshModel( model , new java.io.File( name ) )
			
		}

		/**
		 * Method return highly flexible Gaussian Process to be use in the initial steps with ICP
		 */
		def createGaussianProcess_DefaultKernel( modelReference : ⇒ MeshLandmark ) : 
																						 LowRankGaussianProcess[ _3D, EuclideanVector[ _3D ] ] = {
			
				implicit val rng = scalismo.utils.Random( 42 );

				val zeroMean = scalismo.common.VectorField( RealSpace[_3D] , 
						 										(pt : scalismo.geometry.Point[_3D]) => EuclideanVector.zeros[_3D])

		    val diagKernel1 		= DiagonalKernel[ _3D ]( GaussianKernel( 3 ) * 0.1 , 3 )
		    val diagKernel2 		= DiagonalKernel[ _3D ]( GaussianKernel( 7 ) * 2.0 , 3 )
		    val sum_of_kernels 	= diagKernel1 + diagKernel2
		
		    /*
		     * build a continuous Gaussian Process through GaussianProcess constructor
		     */
		    val gp = GaussianProcess( zeroMean , diagKernel1 )
		
		    //for the reference mesh we will take the firs mesh in the Aligned_dataset
		    val sampler = RandomMeshSampler3D( modelReference.mesh , numberOfPoints = 100 , seed = 45 )

		    /*
		     * In order to be able to sample deformations on a higher resolution mesh, we need a finite 
		     * rank representation of the Gaussian Process
		     * To obtain such a representation in Scalismo, we can use the approximateGP method of the 
		     * LowRankGaussianProcess object
		     */
		    
		    /*
		     * approximateGP has been renamed into approximateGPNystrom
		     * finite rank decomposition of the GP, retaining only the 100 most prominent basis 
		     * functions  
		     */
		    val lowRankGP = scalismo.statisticalmodel.LowRankGaussianProcess.
		    													approximateGPNystrom( gp , sampler , numBasisFunctions = 100 );
				
				return lowRankGP;

		}

		/**
		 * Method return highly flexible Gaussian Process to be use in the initial steps with ICP
		 */
		def createGaussianProcess_HighFlexKernel( modelReference : ⇒ MeshLandmark , scale : Double , 
																							variance : Double , additive : Boolean ) : 
																					LowRankGaussianProcess[ _3D, EuclideanVector[ _3D ] ] = {
			
				// increase flexibility of model to increase fitting to more dynamic mesh
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)

				// defining zero mean
				val zeroMean = scalismo.common.VectorField(RealSpace[_3D], (pt : scalismo.geometry.Point[_3D]) => EuclideanVector.zeros[_3D])

				
				// creates new kernel function with given scale and variance
				val scalarValuedKernel = GaussianKernel[ _3D ]( variance ) * scale;

				// defines symmetry kernel using mirroring at YZ plane
				case class XmirroredKernel( kernel : PDKernel[ _3D ] ) extends PDKernel[ _3D ] {

						override def domain = RealSpace[ _3D ];
					  override def k( x : Point[ _3D ] , y : Point[ _3D ] ) = 
					  																	kernel( Point( x( 0 ) * -1f , x( 1 ) , x( 2 ) ) , y );
				
				}

				// combines kernel with default one to add symmetry correspondence
				// ( instead of mirroring data )
				def symmetrizeKernel( kernel : PDKernel[ _3D ] ) : MatrixValuedPDKernel[ _3D ] = {
				
						val xmirrored = XmirroredKernel( kernel );
					  val k1 		 = DiagonalKernel( kernel , 3 );
					  val k2 		 = DiagonalKernel( xmirrored * -1f , xmirrored , xmirrored );
					  var kcombi = k1 + k2;;
					  if( !additive ){ kcombi = k2; }
					  kcombi;
				
				}

				// creating samples for low rank approximation
				val sampler = RandomMeshSampler3D( modelReference.mesh , numberOfPoints = 300 , seed = 45 )

				// create GP by combine version of the two kernels
				val gp = GaussianProcess[ _3D , EuclideanVector[ _3D ] ](
																								zeroMean , symmetrizeKernel( scalarValuedKernel ) );
				// create low rank approximation from GP using Nystorm approximation
				val lowrankGP = LowRankGaussianProcess.approximateGPNystrom( gp , sampler , 
																																		 numBasisFunctions = 50 )

				// augment given model by (low rank) GP
				return lowrankGP;

		}
		
				/**
		 * Method return highly flexible Gaussian Process to be use in the initial steps with ICP
		 */
		def createGaussianProcess_HighFlexKernelModel( modelReference : ⇒ MeshLandmark , scale : Double , 
																									 variance : Double , additive : Boolean ) : 
																																						StatisticalMeshModel = {
			
				val lowrankGP = createGaussianProcess_HighFlexKernel( modelReference , scale , variance , 
																															additive );

				// augment given model by (low rank) GP
				return StatisticalMeshModel( modelReference.mesh , lowrankGP );

		}

		
		/**
		 * modification kernel to augment a given model
		 * symmetry kernel : 
		 * 		- adds  (  ( X , -Y , -Z ) ( -X , Y , Z) ( -X , Y , Z )  ) symmetry correspondence
		 * 		- increases scalability and flexibility by given scale and variance
		 */
		def modifyModel_SymmetryKernel( model : ⇒ StatisticalMeshModel , scale : Double , 
																		variance : Double , additive : Boolean ) : StatisticalMeshModel = {
			
				// increase flexibility of model to increase fitting to more dynamic mesh
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)
			
				// creates new kernel function with given scale and variance
				val scalarValuedKernel = GaussianKernel[ _3D ]( variance ) * scale;

				// defines symmetry kernel using mirroring at YZ plane
				case class XmirroredKernel( kernel : PDKernel[ _3D ] ) extends PDKernel[ _3D ] {

						override def domain = RealSpace[ _3D ];
					  override def k( x : Point[ _3D ] , y : Point[ _3D ] ) = 
					  																	kernel( Point( x( 0 ) * -1f , x( 1 ) , x( 2 ) ) , y );
				
				}

				// combines kernel with default one to add symmetry correspondence
				// ( instead of mirroring data )
				def symmetrizeKernel( kernel : PDKernel[ _3D ] ) : MatrixValuedPDKernel[ _3D ] = {
				
						val xmirrored = XmirroredKernel( kernel );
					  val k1 		 = DiagonalKernel( kernel , 3 );
					  val k2 		 = DiagonalKernel( xmirrored * -1f , xmirrored , xmirrored );
					  var kcombi = k1 + k2;;
					  if( !additive ){ kcombi = k2; }
					  kcombi;
				
				}

				// create GP by combine version of the two kernels
				val gp = GaussianProcess[ _3D , EuclideanVector[ _3D ] ]( 
																													symmetrizeKernel( scalarValuedKernel ) );
				// create low rank approximation from GP using cholesky decomposition
				val lowrankGP = LowRankGaussianProcess.approximateGPCholesky( 
																											model.referenceMesh.pointSet,
																											gp, relativeTolerance = 0.01, 
						  																				interpolator = NearestNeighborInterpolator() )

				// augment given model by (low rank) GP
				return StatisticalMeshModel.augmentModel( model , lowrankGP );

		}

		/**
		 * modification kernel to augment a given model
		 * scaling kernel : 
		 * 		- homogeneously increase scalability and flexibility of a model by given variance and scale
		 */
		def modifyModel_relaxKernel( model : ⇒ StatisticalMeshModel , scale : Double , 
																 variance : Double ) : StatisticalMeshModel = {
			
				// increase flexibility of model to increase fitting to more dynamic mesh
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random( 42 )
			
				// creates new kernel with given scale and variance
				val relaxKernel = GaussianKernel[ _3D ]( variance ) * scale;
				val kernel 			= DiagonalKernel( relaxKernel , relaxKernel , relaxKernel )
				
				// create GP by combine version of the two kernels
				val gp = GaussianProcess[ _3D , EuclideanVector[ _3D ] ]( kernel );
				// create low rank approximation from GP using cholesky decomposition
				val lowrankGP = LowRankGaussianProcess.approximateGPCholesky( 
																											model.referenceMesh.pointSet ,
																											gp , relativeTolerance = 0.01 , 
						  																				interpolator = NearestNeighborInterpolator() )

				// augment given model by (low rank) GP
				return StatisticalMeshModel.augmentModel( model , lowrankGP );
			
		}

		/**
		 * create model from highly flexible kernel to be use in the initial steps with ICP
		 */
		def buildGaussProcessModel_implicitPosteriori( modelReference : ⇒ MeshLandmark ) 
																																					: StatisticalMeshModel = {

				// create GP from highly flexible kernel																		
		    val lowRankGP = createGaussianProcess_HighFlexKernel( modelReference , 30 , 35 , false );

		    // create model from GP and given reference mesh
				val model = StatisticalMeshModel( modelReference.mesh , lowRankGP );				
								
				return model;
				
		}
		
		/**
		 * Posteriori Model generation using landmarks
		 * The subroutine use the high flexible kernel as a default Prior and train it using posterior
		 * regression. Hereby, only given set of landmark points are considered as correspondence
		 * reference. 
		 */
		def buildGaussProcessModel_Posterior( data : ⇒ MeshLandmarkArrays , 
																					modelReference : ⇒ MeshLandmark ) 
																																					: StatisticalMeshModel = {

				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random( 42 );

				var lowRankGP = createGaussianProcess_HighFlexKernel( modelReference , 30 , 35 , false );
		
				val noise    = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																											 DenseMatrix.eye[ Double ]( 3 )  );

				val refLMPoints = modelReference.lm.map( lm ⇒  lm.point );

				( 0 until data.lms.size ).map( sample ⇒ {

						val lmPoints = data.lms( sample ).map( lm ⇒  lm.point );

						val lmDeform   = ( 0 until lmPoints.size ).map( elem ⇒ { 
							
								( refLMPoints( elem ).toVector - lmPoints( elem ).toVector );							
							
						});
						
						val regressData = ( 0 until lmPoints.size ).map( elem ⇒
			    												IndexedSeq( ( lmPoints( elem ) , lmDeform( elem ) , noise) ) );
					
						regressData.map( trainingElem ⇒ { lowRankGP = lowRankGP.posterior( trainingElem ) } );

				});
				
				val model = StatisticalMeshModel( modelReference.mesh , lowRankGP );				
				
				return model;
				
		}
		
		/**
		 * Posteriori Model generation using uniform sample generator and ICP to find corresponding
		 * points. 
		 *    - using 5000 sample points as correspondence reference
		 *    - as general prior, we use the highly flexible kernel
		 *    - the ICP iteration uses the same model that is going to be trained
		 *      allowing for more flexibility at the beginning and ending with more characteristics
		 *      matching towards the end
		 *    - using a fixed number of ICP iterations
		 */
		def buildGaussProcessModel_PosteriorICP( data : ⇒ MeshLandmarkArrays , 
																						 modelReference : ⇒ MeshLandmark , iterSteps : Int ) 
																																					: StatisticalMeshModel = {

				// iterative ICP non-rigid model fitting using Gaussian Process Regression to train a model
				// from samples using the Posteriori for projection
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)

				// creating reference points to compare the samples with one another 
				val sampler = UniformMeshSampler3D( modelReference.mesh , numberOfPoints = 5000 )

				// extracting pointset from samples
				val points : Seq[ Point[ _3D ] ] = sampler.sample.map( 
																										pointWithProbability ⇒ pointWithProbability._1 );
				// getting the ID of the reference points in the reference mesh																										
				val pointsIDs = points.map( pt ⇒ modelReference.mesh.pointSet.findClosestPoint( pt ).id );	

				//new noise distribution adding additional flexibility to fitting
				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																										DenseMatrix.eye[ Double ]( 3 ) * 0.5 );
							
				println( "sampler primed on reference mesh" )

				var counter : Int = 0;
				
				
			  // creating model using highly flexibel Gaussian Kernel as generalised Prior
				val gp = createGaussianProcess_HighFlexKernel( modelReference , 30 , 35 , false );
				var model = StatisticalMeshModel( modelReference.mesh , gp );
				data.meshes.map( sampleMesh ⇒ { // training looping over the training data sets

		    		println( "performing non-iterative PC for partial shape " + counter )

		    		// creating an overwritable model
						var fittedModel = model;
						var fittedMesh  = model.mean;
		    		var correspondingPoints : Seq[ PointWithId[ _3D ] ] = null;
						
						for( iterstep ← 1 until iterSteps ){

								println( "beginning iteration " + iterstep )

								// generate candidates
								correspondingPoints = pointsIDs.map{ ptID : PointId ⇒ 
									
										val pt 				= fittedMesh.pointSet.point( ptID );
										val ptClosest = sampleMesh.pointSet.findClosestPoint( pt ).point;
										PointWithId( ptClosest , ptID );
									
								}
								println( "found corresponding points to sampler data" + iterstep )
								
								// create successional posterior model
								val iterationTripple = correspondingPoints.map( fittingPoint ⇒ 
																								( fittingPoint.id , fittingPoint.point , noise ) );
								
								fittedModel = model.posterior( iterationTripple.toIndexedSeq );

								println( "regressing to posteriori model done" + iterstep )

								// use found mean as new reference
								fittedMesh = fittedModel.mean;
								/*
								val avgDistance = scalismo.mesh.MeshMetrics.avgDistance( modelReference.mesh, 
																																				 fittedMesh );
								val hausdorffDistance = 
												scalismo.mesh.MeshMetrics.hausdorffDistance( modelReference.mesh, 
																																		 fittedMesh );
		
								println( "Number of iteration: " + iter + ", for target mesh: " + targetNo );
								println( "Average distance: " + avgDistance );
								println( "Hausdorff distance: " + hausdorffDistance + "\n" );
								*/
							
						}
						
						// storing back the posterior from the last found best matching correspondence
						model = fittedModel;
											
				} );
				
				return model;
				
		}

		/**
		 * PCA Model generation using the entire mesh point set as reference assuming global point
		 * ordering in the samples ( applicable using interpolated training sample sets )
		 * PCA uses 10^-6 as relative tolerance
		 */
		def buildGaussProcessModel_PCA( data : ⇒ MeshLandmarkArrays , modelReference : ⇒ MeshLandmark ) 
																																					: StatisticalMeshModel = {

				implicit val rng = scalismo.utils.Random( 42 );

				// create deformation representation for entire training data set
		    val deformation = data.meshes.map( mesh ⇒ {
		    	
						val meshDeform = modelReference.mesh.pointSet.pointIds.map{ id ⇒
							
								mesh.pointSet.point( id ) - modelReference.mesh.pointSet.point( id )														
									    		
						}.toIndexedSeq;
						
						DiscreteField[ _3D , UnstructuredPointsDomain[ _3D ] , EuclideanVector[ _3D ] ](
																												modelReference.mesh.pointSet , meshDeform );
		    	
		    });

		    // interpolate to create continuous deformation field
		    val interp = NearestNeighborInterpolator[ _3D , EuclideanVector[ _3D ] ]()

		    val contDeform = deformation.map( meshDeform ⇒ meshDeform.interpolate( interp ) );

		    // use PCA on deformation set using the given reference
		    val lowRankGaussian = DiscreteLowRankGaussianProcess.createUsingPCA( 
																	    										modelReference.mesh.pointSet , contDeform , 
																	    										RelativeTolerance( 1e-6 ) );

		    // create model using same reference the the (low rank) GP from PCA
				val model = StatisticalMeshModel( modelReference.mesh , lowRankGaussian );
				
				return model;
				
		}

		/**
		 * PCA Model generation using uniform sample generator and ICP to find corresponding points. 
		 *    - using 5000 sample points as correspondence reference
		 *    - for ICP we use the model from the the highly flexible kernel
		 *    - using a fixed number of ICP iterations
		 *    - the found corresponding mesh points form the ICP are used as reference points 
		 *    	for the PCA
		 *    - PCA uses 10^-6 as relative tolerance
		 */
		def buildGaussProcessModel_PCA_IPC( data : ⇒ MeshLandmarkArrays , 
																				modelReference : ⇒ MeshLandmark , iterSteps : Int ) : 
																																						StatisticalMeshModel = {

				// iterative ICP non-rigid model fitting using Gaussian Process Regression to train a model
				// from samples using PCA
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)
				
				// creating reference point set from uniform mesh samples
				var points : Seq[ Point[ _3D ] ] = null;
				{
						val sampler = UniformMeshSampler3D( modelReference.mesh , numberOfPoints = 5000 )

						points = sampler.sample.map( pointWithProbability ⇒ pointWithProbability._1 );
				}
				// find corresponding point IDs in reference mesh 
				val pointsIDs = points.map( pt ⇒ modelReference.mesh.pointSet.findClosestPoint( pt ).id );	

				//new noise distribution adding additional flexibility to fitting
				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																														DenseMatrix.eye[ Double ]( 3 ) * 0.5 );
							
				println( "reference mesh points generated" )

				var modelCorrespondingPoints : Seq[ Seq[ PointWithId[ _3D] ] ] = null;
				{ // localise IPC variables
					var counter : Int = 0;

			    	// creating model using highly flexibel Gaussian Kernel as generalised Prior
					val gp = createGaussianProcess_HighFlexKernel( modelReference , 30 , 35 , false );
					var model = StatisticalMeshModel( modelReference.mesh , gp );
					// looping over the training data sets to find corresponding points for all samples
					modelCorrespondingPoints = data.meshes.map( sampleMesh ⇒ {
	
			    		println( "performing non-iterative PC for partial shape " + counter )

			    		// re-initialising mean for each sample
							var fittedMesh  = model.mean;
			    		var correspondingPoints : Seq[ PointWithId[ _3D ] ] = null;
							
							for( iterstep ← 0 until iterSteps ){
	
									println( "beginning iteration " + iterstep )
	
									// generate candidates
									correspondingPoints = pointsIDs.map{ ptID : PointId ⇒ 
										
											val pt 				= fittedMesh.pointSet.point( ptID );
											val ptClosest = sampleMesh.pointSet.findClosestPoint( pt ).point;
											PointWithId( ptClosest , ptID );
										
									}
									println( "found corresponding points to sampler data" + iterstep )
									
									// create successional posterior model
									val iterationTripple = correspondingPoints.map( fittingPoint ⇒ 
																									( fittingPoint.id , fittingPoint.point , noise ) );
									
									val fittedModel = model.posterior( iterationTripple.toIndexedSeq );
	
									println( "regressing to posteriori model done" + iterstep )

									// use found mean as new reference
									fittedMesh = fittedModel.mean;
									/*
									val avgDistance = scalismo.mesh.MeshMetrics.avgDistance( modelReference.mesh, 
																																					 fittedMesh );
									val hausdorffDistance = 
													scalismo.mesh.MeshMetrics.hausdorffDistance( modelReference.mesh, 
																																			 fittedMesh );
			
									println( "Number of iteration: " + iter + ", for target mesh: " + targetNo );
									println( "Average distance: " + avgDistance );
									println( "Hausdorff distance: " + hausdorffDistance + "\n" );
									*/
								
							}
	
							// storing back the last found ( and best matching ) correspondence for each sample
							correspondingPoints;
												
					} );
				}
				
				// create deformation fields for all samples using the reference points 
				// with the found correspondences sample-wise
		    val deformation = ( 0 until modelCorrespondingPoints.size ).map( meshNo ⇒ {
		    	
						val meshDeform = ( 0  until modelReference.mesh.pointSet.numberOfPoints ).map{ elem ⇒
							
								( modelCorrespondingPoints( meshNo )( elem ).point - points( elem ) )
									    		
						}.toIndexedSeq;
						
						DiscreteField[ _3D , UnstructuredPointsDomain[ _3D ] , EuclideanVector[ _3D ] ](
																												modelReference.mesh.pointSet , meshDeform );
		    	
		    });

		    // interpolate to create continuous deformation field
		    val interp = NearestNeighborInterpolator[ _3D , EuclideanVector[ _3D ] ]()

		    val contDeform = deformation.map( meshDeform ⇒ meshDeform.interpolate( interp ) );
		    
		    // use PCA on deformation set using the original reference
		    val lowRankGaussian = DiscreteLowRankGaussianProcess.createUsingPCA( 
		    										modelReference.mesh.pointSet , contDeform , RelativeTolerance( 1e-6 ) );
		    
		    // create model using same reference the the (low rank) GP from PCA
				val model = StatisticalMeshModel( modelReference.mesh , lowRankGaussian );
				
				return model;
				
		};

		/**
		 * PCA Model generation using uniform sample generator and ICP to find interpolated sample 
		 * representation and corresponding points. 
		 *    - using 5000 sample points as correspondence reference
		 *    - for ICP we use the model from the the highly flexible kernel
		 *    - using a fixed number of ICP iterations
		 *    - the found posterior representation of each sample will be used as interpolated sample
		 *    	set ( considering the most likely candidate, i.e. mean )
		 *    - PCA uses 10^-6 as relative tolerance
		 */
		def buildGaussProcessModel_ICP_relaxed_PCA( data : ⇒ MeshLandmarkArrays , 
																			  			  modelReference : ⇒ MeshLandmark , iterSteps : Int ) : 
																			  			  														StatisticalMeshModel = {

				// iterative ICP non-rigid model fitting using Gaussian Process Regression to train a model
				// from samples using PCA
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)
				
				// creating reference points to compare the samples with one another 
				var points : Seq[ Point[ _3D ] ] = null;
				{
						val sampler = UniformMeshSampler3D( modelReference.mesh , numberOfPoints = 5000 )

						points = sampler.sample.map( pointWithProbability ⇒ pointWithProbability._1 );
				}
				// getting the ID of the reference points in the reference mesh																										
				val pointsIDs = points.map( pt ⇒ modelReference.mesh.pointSet.findClosestPoint( pt ).id );	

				//new noise distribution adding additional flexibility to fitting
				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																														DenseMatrix.eye[ Double ]( 3 ) * 0.5 );
							
				println( "sampler primed on reference mesh" )

				var ICP_relaxed : Seq[ TriangleMesh[ _3D ] ] = null;
				{ // localise IPC variables
						var counter : Int = 0;
						
			    	// creating model using highly flexibel Gaussian Kernel as generalised Prior
						val gp 		= createGaussianProcess_HighFlexKernel( modelReference , 30 , 35 , false );
						var model = StatisticalMeshModel( modelReference.mesh , gp );
						// training looping over the training data sets
						ICP_relaxed = data.meshes.map( sampleMesh ⇒ { 
		
				    		println( "performing non-iterative PC for partial shape " + counter )
				    		
		    				// creating an overwritable model
								var fittedMesh  = model.mean;
				    		var correspondingPoints : Seq[ PointWithId[ _3D ] ] = null;
								
								for( iterstep ← 0 until iterSteps ){
		
										println( "beginning iteration " + iterstep )
		
										// generate candidates
										correspondingPoints = pointsIDs.map{ ptID : PointId ⇒ 
											
												val pt 				= fittedMesh.pointSet.point( ptID );
												val ptClosest = sampleMesh.pointSet.findClosestPoint( pt ).point;
												PointWithId( ptClosest , ptID );
											
										}
										println( "found corresponding points to sampler data" + iterstep )
										
										// create successional posterior model
										val iterationTripple = correspondingPoints.map( fittingPoint ⇒ 
																										( fittingPoint.id , fittingPoint.point , noise ) );
										
										val fittedModel = model.posterior( iterationTripple.toIndexedSeq );
		
										println( "regressing to posteriori model done" + iterstep )
	
										// use found mean as new reference
										fittedMesh = fittedModel.mean;
										/*
										val avgDistance = scalismo.mesh.MeshMetrics.avgDistance( modelReference.mesh, 
																																						 fittedMesh );
										val hausdorffDistance = 
														scalismo.mesh.MeshMetrics.hausdorffDistance( modelReference.mesh, 
																																				 fittedMesh );
				
										println( "Number of iteration: " + iter + ", for target mesh: " + targetNo );
										println( "Average distance: " + avgDistance );
										println( "Hausdorff distance: " + hausdorffDistance + "\n" );
										*/
									
								}
		
								// storing back the posterior from the last found best matching correspondence
								fittedMesh;
													
						} );
				}
				
				// create deformation fields for all samples using the reference points 
				// with the found correspondences and the interpolated mesh sample-wise
		    val deformation = ( 0 until ICP_relaxed.size ).map( mesh ⇒ {
		    	
						val meshDeform = ( 0 until pointsIDs.size ).map{ elem ⇒
							
								( ICP_relaxed( mesh ).pointSet.point( pointsIDs( elem ) ) - points( elem ) )														
									    		
						}.toIndexedSeq;
						
						DiscreteField[ _3D , UnstructuredPointsDomain[ _3D ] , EuclideanVector[ _3D ] ](
																												modelReference.mesh.pointSet , meshDeform );
		    	
		    });
		    
		    // interpolate to create continuous deformation field
		    val interp = NearestNeighborInterpolator[ _3D , EuclideanVector[ _3D ] ]()

		    val contDeform = deformation.map( meshDeform ⇒ meshDeform.interpolate( interp ) );
		    
		    // use PCA on deformation set using the original reference
		    val lowRankGaussian = DiscreteLowRankGaussianProcess.createUsingPCA( 
		    										modelReference.mesh.pointSet , contDeform , RelativeTolerance( 1e-6 ) );
		    
		    // create model using same reference the the (low rank) GP from PCA
				val model = StatisticalMeshModel( modelReference.mesh , lowRankGaussian );
				
				return model;
				
		}

		/**
		 * PCA Model generation using uniform sample generator and ICP to find interpolated sample 
		 * representation and corresponding points. 
		 *    - using 5000 sample points as correspondence reference
		 *    - for ICP we use the model from the the highly flexible kernel
		 *    - using a fixed number of ICP iterations
		 *    - the found posterior representation of each sample will be used as interpolated sample
		 *    	set ( considering the most likely candidate, i.e. mean )
		 *    - PCA uses 10^-6 as relative tolerance
		 */
		def buildGaussProcessModel_ICP_relaxed_PCA2( data : ⇒ MeshLandmarkArrays , 
																			  			  modelReference : ⇒ MeshLandmark , iterSteps : Int ) : 
																			  			  														StatisticalMeshModel = {

				// iterative ICP non-rigid model fitting using Gaussian Process Regression to train a model
				// from samples using PCA
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)
				
				var prior_model = createGaussianProcess_HighFlexKernelModel( modelReference, 30, 35, false );

				//new noise distribution adding additional flexibility to fitting
				val noiseParameter : Double = 1
				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																										noiseParameter * DenseMatrix.eye[ Double ]( 3 ) );
							
				println( "sampler primed on reference mesh" )

				var ICP_relaxed : Seq[ TriangleMesh[ _3D ] ] = null;
				{ // localise IPC variables
						var counter : Int = 0;

						// training looping over the training data sets
						ICP_relaxed = data.meshes.map( sampleMesh ⇒ { 
		
				    		println( "performing non-iterative PC for partial shape " + counter )
				    		
		    				// creating an overwritable model
								var meanMesh  = modelReference.mesh;
				    		
				    		//ICP iteration
								for( iterstep ← 0 until iterSteps ){
		
										println( "beginning iteration " + iterstep )
		
										// generate candidates
										val correspondingPoints = meanMesh.pointSet.pointIds.toIndexedSeq.map{ 
																																									ptID : PointId ⇒ 
											
												val pt 				= meanMesh.pointSet.point( ptID );
												val ptClosest = sampleMesh.pointSet.findClosestPoint( pt ).point;
												PointWithId( ptClosest , ptID );
											
										} 
										println( "found corresponding points to sampler data" + iterstep )
										
										// create successional posterior model
										val iterationTripple = correspondingPoints.map( fittingPoint ⇒ 
																										( fittingPoint.id , fittingPoint.point , noise ) );
										
										val fittedModel = prior_model.posterior( iterationTripple.toIndexedSeq );
		
										println( "regressing to posteriori model done" + iterstep )
	
										// use found mean as new reference
										meanMesh = fittedModel.mean;
										/*
										val avgDistance = scalismo.mesh.MeshMetrics.avgDistance( modelReference.mesh, 
																																						 fittedMesh );
										val hausdorffDistance = 
														scalismo.mesh.MeshMetrics.hausdorffDistance( modelReference.mesh, 
																																				 fittedMesh );
				
										println( "Number of iteration: " + iter + ", for target mesh: " + targetNo );
										println( "Average distance: " + avgDistance );
										println( "Hausdorff distance: " + hausdorffDistance + "\n" );
										*/
									
								}
		
								// storing back the posterior from the last found best matching correspondence
								meanMesh;
													
						} );
				}
				
				// create deformation fields for all samples using the reference points 
				// with the found correspondences and the interpolated mesh sample-wise
		    val deformation = ( 0 until ICP_relaxed.size ).map( mesh ⇒ {
		    	
						val meshDeform = modelReference.mesh.pointSet.pointIds.toIndexedSeq.map{ id ⇒
							
								( ICP_relaxed( mesh ).pointSet.point( id ) - 
																													modelReference.mesh.pointSet.point( id ) )														
									    		
						}.toIndexedSeq;
						
						DiscreteField[ _3D , UnstructuredPointsDomain[ _3D ] , EuclideanVector[ _3D ] ](
																												modelReference.mesh.pointSet , meshDeform );
		    	
		    });
		    
		    // interpolate to create continuous deformation field
		    val interp = NearestNeighborInterpolator[ _3D , EuclideanVector[ _3D ] ]()

		    val contDeform = deformation.map( meshDeform ⇒ meshDeform.interpolate( interp ) );
		    
		    // use PCA on deformation set using the original reference
		    val lowRankGaussian = DiscreteLowRankGaussianProcess.createUsingPCA( 
		    										modelReference.mesh.pointSet , contDeform , RelativeTolerance( 1e-6 ) );
		    
		    // create model using same reference the the (low rank) GP from PCA
				val model = StatisticalMeshModel( modelReference.mesh , lowRankGaussian );
				
				return model;
				
		}
		
		/**
		 * Returns whether given point is at the border of a given mesh. 
		 * Here, we assume TriangleMesh with each point having upto 6 neighbours (not a border)
		 */
		def isMeshBorderPoint( mesh : ⇒ TriangleMesh[ _3D ] , pt : PointId ) : Boolean = {
			
				val connections = mesh.triangulation.adjacentTrianglesForPoint( pt );
				
				/* 
				 * in a triangular mesh, each mesh point that is not part of a border ( even corner ) has
				 * 6 neighbours
				 */
				if( connections.length < 6 ){
					return true;
				} else {
					return false;
				}
			
		}

		/**
		 * Find maximal local the spacing distance inbetween mesh points at a given point (pt).
		 */
		def maxMeshSubspace2( pt : PointId , mesh : ⇒ TriangleMesh[ _3D ] ) : Double = {

				val origin		 = mesh.pointSet.point( pt );
				val neighbours = mesh.triangulation.adjacentPointsForPoint( pt )
				var dist : Double = 0;
				for( neighbour ← neighbours ){ // loop over all neighbours to find local maximal grid space
					
						val newdist = ( mesh.pointSet.point( neighbour ) - origin ).norm2;
						if( dist < newdist ) dist = newdist;
					
				}
				return dist;
			
		}

		/**
		 * Find maximal local the spacing distance inbetween mesh points at a given point (pt) 
		 * considering only mesh-edges that are at the border of the mesh. 
		 * Here, we assume TriangleMesh with each point having upto 6 neighbours (not a border)
		 */
		def maxMeshSubspace2Border( pt : PointId , mesh : ⇒ TriangleMesh[ _3D ] ) : Double = {

				val origin		 = mesh.pointSet.point( pt );
				val neighbours = mesh.triangulation.adjacentPointsForPoint( pt )
				var dist : Double = 0;
				for( neighbour ← neighbours ){ // loop over neighbours
					
						val newdist = ( mesh.pointSet.point( neighbour ) - origin ).norm2;
						val connections = mesh.triangulation.adjacentPointsForPoint( neighbour );
						//only consider grid spacing along border
						if( connections.length < 6 && dist < newdist ) dist = newdist;
					
				}
				return dist;
			
		}
		
		/**
		 * Method returns whether or not the new point of the new border mesh actually defines a new
		 * border or is overlapping with the mesh defining the old border. 
		 * ⇒  this is used to evaluate the extend of the missing information in a mesh with gaps ( mesh
		 * defining the old border ) and is indeed present int the complete mesh ( mesh defining the new
		 * border )
		 * The two points are considered to be each the closest one to the other within the other mesh. 
		 * ( given as argument to reduce unnecessary recalculation of the closest points )
		 * Here, we assume TriangleMesh with each point having upto 6 neighbours (not a border)
		 */
		def isDefinesNewBorder( pt_new_closest : PointId , meshNewBorder : ⇒ TriangleMesh[ _3D ] , 
														pt_old_closest : PointId , meshOldBorder : ⇒ TriangleMesh[ _3D ] ) : 
																																												Boolean = {
			
				// get neighbour points
				val neighbours	= meshNewBorder.triangulation.adjacentPointsForPoint( pt_new_closest );
				val point1			= meshNewBorder.pointSet.point( pt_new_closest );
				
				// if distance between corrsponding points is smaller grid space
				var smaller : Boolean 	 = true;	
				/*
				 *  if one neighbour lies in inside the mesh of the completed missing information
				 *  ( compactness criteria of missing information pieces (gap completions) )
				 */
				var innerPoint : Boolean = false;

				val closestOld1	= meshOldBorder.pointSet.findClosestPoint( point1 );
				val dist1       = ( point1 - closestOld1.point ).norm2;
				for( neighbour ← neighbours ){ // loop over neighbours
					
						val point2			= meshNewBorder.pointSet.point( neighbour );
						val distNew			= ( point1 - point2 ).norm2; 														// partial mesh space
//						val closestOld2	= meshOldBorder.pointSet.findClosestPoint( point2 );	
//						var distOld			= ( closestOld1.point - closestOld2.point ).norm2;		// complete mesh space
						if( dist1 < distNew ) { smaller = false; } // has farther mesh point than corresp dist
						val point2Cor   = meshOldBorder.pointSet.findClosestPoint( point2 );
						if( !isMeshBorderPoint( meshNewBorder , neighbour ) &&  // has inside neighbour
								isMeshBorderPoint( meshOldBorder , point2Cor.id ) ){ innerPoint = true; }
						
				}
				if( smaller ) return true;		// all mesh spaces are smaller than distance to border
				
				if( innerPoint ) return true;// is directly connected to guaranteed member ⇒ also member
				
				return false;			
			
		}
		
		/**
		 * Method finds corresponding elements on either meshmodel and reduces the complete model
		 * the the missing mesh information of the model with gaps
		 */
		def reduceToGap( meshWithGap : ⇒ TriangleMesh[ _3D ] , meshComplete : ⇒ TriangleMesh[ _3D ] ) :
																																						Seq[ PointId ] = {

				// get list of borderpoints of the mesh with gaps where the complete mesh has no border
				var gapFillingPointIDs : Seq[ PointId ] = ( 0 until meshComplete.pointSet.numberOfPoints ).
					map( k ⇒  { // loop over all mesh points as potential candidates

							val pointMeshComplete = meshComplete.pointSet.points.toIndexedSeq( k );
							val pointMeshGap 			= meshWithGap.pointSet.findClosestPoint( pointMeshComplete );
							if( isMeshBorderPoint( meshWithGap , pointMeshGap.id ) ){
								
									if( !isMeshBorderPoint( meshComplete , PointId( k ) ) ){
										
											pointMeshGap.id; // fulfilled first criterion
										
									} else if ( isDefinesNewBorder( PointId( k ) , meshComplete , 
																									pointMeshGap.id , meshWithGap ) ) {
										
											pointMeshGap.id; // is valid border case : fulfilled 2nd or 3rd criterion
										
									} else { PointId( -1 ); }
								
							} else { PointId( -1 ); }
							
					} );
				// filter to valid ID (verified targets) 
				gapFillingPointIDs 		= gapFillingPointIDs.filter( _.id >= 0 ); 
												
				return gapFillingPointIDs;
							
		}

		/**
		 * projection of model to partial meshes by (partial) posterior regression using landmarks
		 * non-iterative CP version assuming global point ordering in the samples 
		 * ( applicable using landmarks or interpolated training sample sets )
		 */
		def projectPartialMesh_nonIterative_CP(  model : ⇒ StatisticalMeshModel , 
																						 data : ⇒ MeshLandmarkArrays , 
																						 modelReference : ⇒ MeshLandmark ) : 
																																	Array[ TriangleMesh[ _3D ] ] = {
			
				// non Iterative CP using Gaussian Process Regression to complete partial shape
				// ... from Scalismo Tutorials
				val newModel = modifyModel_relaxKernel( model , 10.0 , 70.0 );
				
				var counter : Int = 0;

				val noise    = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																											 DenseMatrix.eye[ Double ]( 3 ) * 0.5 );

				//loop over all partial meshes
				val bestFitMissingPiece : IndexedSeq[ TriangleMesh[ _3D ] ] = 
																											( 0 until data.meshes.size ).map( elem ⇒ {
					
		    		println( "performing non-iterative PC for partial shape " + counter )

		    		// use landmarks as reference points
		    		val modelRefLMPoints : Seq[ Point[ _3D ] ]	= modelReference.lm.map( lm ⇒ lm.point );
						val LMsPoints 			 : Seq[ Point[ _3D ] ]	= data.lms( elem ).map( lm ⇒ lm.point );
						
						// create deformation field representation
						val domain		= UnstructuredPointsDomain( LMsPoints.toIndexedSeq );
						val deform		= ( 0 until modelRefLMPoints.size ).
																										map( p ⇒ LMsPoints( p ) - modelRefLMPoints( p ) );
						val	defField  = DiscreteField[ _3D , UnstructuredPointsDomain[ _3D ] , 
																					 EuclideanVector[ _3D ] ]( domain , deform );

						println( "found deformation field " )

						// perform posteriori regression
						val Regress	  = for( ( ref , target ) ← modelRefLMPoints zip LMsPoints ) yield { 
							
								val refPointID = newModel.referenceMesh.pointSet.findClosestPoint( ref ).id;
								( refPointID , target , noise )
			
						}
						val posti : StatisticalMeshModel = newModel.posterior( Regress.toIndexedSeq );

						println( "regression to posteriori model done" )
						
						// restrict model to target-area
						val gapPointIDs = reduceToGap( data.meshes( elem ) , posti.mean );

						println( "missing parts in partial shapes identified" )

						// trim model to target point set
						val postiRestr : StatisticalMeshModel = posti.marginal( gapPointIDs.toIndexedSeq );
						
						//get most likely sample
						val gapFilling : TriangleMesh[ _3D ] = postiRestr.mean;

						println( "model trimed and most likely shape evaluated" )
						
						counter = counter + 1;
				
						//return most likely sample as result
						gapFilling
						
				} );
				
				return bestFitMissingPiece.toArray;
				    }
		
		/**
		 * projection of model to partial meshes by (partial) posterior regression using IPC
		 *    - using 5000 sample points as correspondence reference
		 *    - using a fixed number of ICP iterations
		 */
		def projectPartialMesh_Iterative_CP(  model : ⇒ StatisticalMeshModel , 
																				  data : ⇒ MeshLandmarkArrays , 
																					modelReference : ⇒ MeshLandmark , 
																					numberOfIterations : Int ) : 
																																	Array[ TriangleMesh[ _3D ] ] = {
			
				// iterative ICP non-rigid model fitting using Gaussian Process Regression to complete 
				//	Partial Shapes
				// ... from Scalismo Tutorials
				implicit val rng = scalismo.utils.Random(42)
				
				// creating reference point set from uniform mesh samples
				var points : Seq[ Point[ _3D ] ] = null;
				{
						val sampler = UniformMeshSampler3D( modelReference.mesh , numberOfPoints = 5000 )

						points = sampler.sample.map( pointWithProbability ⇒ pointWithProbability._1 );
				}
				// find corresponding point IDs in reference mesh 
				val pointsIDs = points.map( pt ⇒ modelReference.mesh.pointSet.findClosestPoint( pt ).id );	

				println( "reference mesh points generated" )

				var counter : Int = 0;

				//noise distribution adding additional flexibility to fitting
				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ) , 
																										DenseMatrix.eye[ Double ]( 3 ) * 0.5 );
							

				// looping over the training data sets to find corresponding points for all samples
				val fittedMeshes = data.meshes.map( targetmesh ⇒ {

		    		println( "performing non-iterative PC for partial shape " + counter )
		    		
			    	// creating an overwritable model
						var fittedMesh  = model.mean;
						var fittedModel = model;
					
						for( iterstep ← 1 until numberOfIterations ){

								println( "beginning iteration " + iterstep )

								// generate candidates
								val correspondingPoints : Seq[ PointWithId[ _3D ] ] = pointsIDs.map{ ptID : PointId ⇒ 
									
										val pt 				= fittedMesh.pointSet.point( ptID );
										val ptClosest = targetmesh.pointSet.findClosestPoint( pt ).point;
										PointWithId( ptClosest , ptID );
									
								}
								println( "found corresponding points to sampler data" + iterstep )
								
								// create successional posterior model
								val iterationTripple = correspondingPoints.map( fittingPoint ⇒ 
																								( fittingPoint.id , fittingPoint.point , noise ) );
								
								fittedModel = model.posterior( iterationTripple.toIndexedSeq );

								println( "regressing to posteriori model done" + iterstep )

								// use found mean as new reference
								fittedMesh = fittedModel.mean;
							
						}

						// restrict model to target-area
						val gapPointIDs = reduceToGap( targetmesh , fittedMesh );
						
						// trim model to target point set
						val postiRestr : StatisticalMeshModel = fittedModel.marginal( gapPointIDs.toIndexedSeq );
						
						println( "model trimed and most likely shape evaluated" )

						counter = counter + 1;
						
						//return most likely sample as result
						fittedModel.mean;
					
				} );
				
				return fittedMeshes;
				
		}
		
		/**
		 * projection of model to partial meshes by (partial) posterior regression using ICP on 
		 * reference landmarks
		 *    - using a fixed number of ICP iterations
		 */
		def projectPartialMesh_RawGP_Iterative_CP( 
																				gp : ⇒ LowRankGaussianProcess[_3D, EuclideanVector[ _3D ] ] , 
																				data : ⇒ MeshLandmarkArrays , 
																				modelReference : ⇒ MeshLandmark , 
																				numberOfIterations : Int , tempLoc : String ) : 
																																		Array[ TriangleMesh[ _3D ] ] = {

				implicit val rng = scalismo.utils.Random(42)

				// generate reference point set using landmarks and find corresponding reference IDs
				val modelReferencePoints : Seq[ Point[ _3D ] ] = modelReference.lm.map( lm ⇒ lm.point )
				val modelReferenceIDs    = modelReferencePoints.map( pt ⇒ 
																						{ modelReference.mesh.pointSet.findClosestPoint( pt ) } );
				
				val modelReferencePointsDomain = UnstructuredPointsDomain( modelReferencePoints.toIndexedSeq );
																		
		
		    //___________________________Using the candidate deformation fields for regression________________________________//
		
		    val littleNoise : MultivariateNormalDistribution = MultivariateNormalDistribution( 
   												DenseVector.zeros[ Double ]( 3 ) , 0.1 * DenseMatrix.eye[ Double ]( 3 ) )
		    println("NOISEEEE")
				
		    println( "START iteration" )

		    val currentPoints = modelReferencePoints;
				
				val nbIterations : Int = 1;
		    
				//loop over all partial meshes
		    val results : IndexedSeq[ TriangleMesh[ _3D ] ] = ( 0 until data.meshes.size ).
																																									map( targetNo ⇒ {
		    	
			    // re-initialise reference for each sample
			    var iter = 0
			    var mean = modelReference.mesh;
			    while( iter < nbIterations ){
			    	
			    		// find candidates of corresponding points
				      val candidates = currentPoints.map{pt => data.meshes( targetNo ).pointSet.
				      																												findClosestPoint( pt ).point }

				      // create deformation representation
							val deff = ( 0 until modelReference.lm.size ).map( j => ( candidates( j ) - 
																																				modelReferencePoints( j ) ) )
																																				
							// create regression data set
				      val trainData : IndexedSeq[ ( Point[ _3D ] , EuclideanVector[ _3D ] , 
				      										 MultivariateNormalDistribution ) ] =  
			      											 				(modelReferencePoints zip deff).map{ case (dp, deform) => 
		      											 																		( dp , deform , littleNoise ) 
		      											 					}.toIndexedSeq;
		      											 					
							// create successional posterior model
				      val posteriorGP = gp.posterior( trainData )
				      val  deformedModel = StatisticalMeshModel(modelReference.mesh, posteriorGP)
				
							//get most likely sample
							mean = deformedModel.mean 
				      
							iter = iter + 1;		    	

		    		  println("Number of iteration: "+ iter + ", for target mesh: " + targetNo)

			    }
		    
					//return most likely sample as result
		    	mean;
		    	
		    })

		   	println(" END iteration")
		   	
		   	return results.toArray;
			
		}
		
		/**
		 * projection of model to partial meshes by (partial) posterior regression using IPC
		 *    - using 5000 sample points as correspondence reference
		 *    - using a fixed number of ICP iterations
		 *    - evaluate ICP result using average distance and Hausdorff distance
		 */
		def projectPartialMesh_Iterative_CP_optim( shapeModel : ⇒ StatisticalMeshModel , 
																					  	 data : ⇒ MeshLandmarkArrays , 
																							 modelReference : ⇒ MeshLandmark , 
																							 numberOfIterations : Int , tempLoc : String ) 
																																	: Array[ TriangleMesh[ _3D ] ] = {

		    // setting a seed for the random generator to allow for reproducible results
		    implicit val rng = scalismo.utils.Random( 42 )
		
		    //___________________________Recursion for Iteration________________________________//

		    println("START RECURSION")
		    
				// creating reference point set from uniform mesh samples
		    var newLandmarks_ref_id : Seq[ PointId ] = null;
		    {
		    		val newLandmarks_ref = UniformMeshSampler3D( modelReference.mesh , 1000 ).sample.
		    																																					map( s => s._1 );
		    		newLandmarks_ref_id	 = newLandmarks_ref.map( pt ⇒ 
		    																		modelReference.mesh.pointSet.findClosestPoint( pt ).id );
		    }

		    var noiseParameter : Double = 0;

				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ),
																										noiseParameter * DenseMatrix.eye[ Double ]( 3 ) );
		    
				val numberOfIter : Int = 5;
				// looping over the training data sets
				val meshesCompleted = ( 0 until data.meshes.size ).map( targetNo ⇒ {

						// scan for pre-evaluated results sample-wise and skip evaluation of such partials
						var posterioriModel 	 : StatisticalMeshModel = null;
						var reference_iterated : TriangleMesh[_3D ] 	= null;
						posterioriModel = readModelFile( tempLoc + "/partialModel_" + targetNo + ".h5" )
						if( posterioriModel == null ){ // only evaluate partials without existing solution
							
			    	// creating an overwritable model
						reference_iterated 	= modelReference.mesh;
						posterioriModel 	 	= shapeModel;
						( 0 until numberOfIter ).map( iter ⇒ {
		
								//compute candidates , returnes ids +  new candidates on the target Mesh!
								val candidates = newLandmarks_ref_id.map( id => {
			
										val pt = reference_iterated.pointSet.point( id );
										val clostestPtsOnTarget = 
																				data.meshes( targetNo ).pointSet.findClosestPoint( pt ).point;
										( id, clostestPtsOnTarget )
			
								} );
			
								val newCandidates = candidates.map( c => c._2 ).toIndexedSeq;
			
								// create successional posterior model
								val trainData = candidates.map( c => ( c._1 , c._2 , noise ) );
			
								posterioriModel = shapeModel.posterior( trainData.toIndexedSeq );
			
								//get most likely sample
								reference_iterated = posterioriModel.mean;
		
								// evaluate correspondence
								val avgDistance = scalismo.mesh.MeshMetrics.avgDistance(	modelReference.mesh, 
																																					reference_iterated );
								val hausdorffDistance = 
													scalismo.mesh.MeshMetrics.hausdorffDistance( modelReference.mesh, 
																																			 reference_iterated );
			
								println( "Number of iteration: " + iter + ", for target mesh: " + targetNo );
								println( "Average distance: " + avgDistance );
								println( "Hausdorff distance: " + hausdorffDistance + "\n" );
		
						} );
					
						writeModelFile( posterioriModel , tempLoc + "/partialModel_" + targetNo + ".h5" )
						}
					
						//return most likely sample as result
						posterioriModel.mean
						
			  } ).toArray;

			  return meshesCompleted;					      

		}
		
		/**
		 * projection of model to partial meshes by (partial) posterior regression using ICP on 
		 * reference landmarks
		 * UI version returning graphical feedback to user
		 *    - using a fixed number of ICP iterations
		 */
		def projectPartialMesh_Iterative_CP_optim_UI( shapeModel : ⇒ StatisticalMeshModel , 
																					  	 		data : ⇒ MeshLandmarkArrays , 
																							 		modelReference : ⇒ MeshLandmark , 
																							 		numberOfIterations : Int , 
																							 		UI : ⇒ ScalismoUI , tempLoc : String ) : 
																							 											Array[ TriangleMesh[ _3D ] ] = {

				// iterative ICP non-rigid model fitting using Gaussian Process Regression to complete 
				//	Partial Shapes
				// ... from Scalismo Tutorials
		    implicit val rng = scalismo.utils.Random( 42 )
		
		    //___________________________Recursion for Iteration________________________________//

		    println("START RECURSION")
		    
				// creating reference point set from uniform mesh samples
		    var newLandmarks_ref_id : Seq[ PointId ] = null;
		    {

		    		val newLandmarks_ref = UniformMeshSampler3D( modelReference.mesh , 1000 ).sample.
		    																																					map( s => s._1 );
		    		newLandmarks_ref_id	 = newLandmarks_ref.map( pt ⇒ 
		    																		modelReference.mesh.pointSet.findClosestPoint( pt ).id );
		    }

				//noise distribution adding additional flexibility to fitting
		    var noiseParameter : Double = 0;
				val noise = MultivariateNormalDistribution( DenseVector.zeros[ Double ]( 3 ),
																										noiseParameter * DenseMatrix.eye[ Double ]( 3 ) );

								
				val numberOfIter : Int = 5;
				// looping over the training data sets to find corresponding points for all samples
				val meshesCompleted = ( 0 until data.meshes.size ).map( targetNo ⇒ {

						// scan for pre-evaluated results sample-wise and skip evaluation of such partials
						var posterioriModel 	 : StatisticalMeshModel = null;
						var reference_iterated : TriangleMesh[_3D ] 	= null;
						posterioriModel = readModelFile( tempLoc + "/partialModel_" + targetNo + ".h5" )
						if( posterioriModel == null ){ // only evaluate partials without existing solution
							
			    	// creating an overwritable model
						reference_iterated 	= modelReference.mesh;
						posterioriModel 	 	= shapeModel;
						( 0 until numberOfIter ).map( iter ⇒ {
		
								//compute candidates , returnes ids +  new candidates on the target Mesh!
								val candidates = newLandmarks_ref_id.map( id => {
			
										val pt = reference_iterated.pointSet.point( id );
										val clostestPtsOnTarget = 
																				data.meshes( targetNo ).pointSet.findClosestPoint( pt ).point;
										( id, clostestPtsOnTarget )
			
								} );
			
								val newCandidates = candidates.map( c => c._2 ).toIndexedSeq;
			
								// create successional posterior model
								val trainData = candidates.map( c => ( c._1 , c._2 , noise ) );
			
								posterioriModel = shapeModel.posterior( trainData.toIndexedSeq );
			
								//get most likely sample
								reference_iterated = posterioriModel.mean;
		
								// evaluate correspondence
								val avgDistance = scalismo.mesh.MeshMetrics.avgDistance(	modelReference.mesh, 
																																					reference_iterated );
								val hausdorffDistance = 
													scalismo.mesh.MeshMetrics.hausdorffDistance( modelReference.mesh, 
																																			 reference_iterated );
			
								println( "Number of iteration: " + iter + ", for target mesh: " + targetNo );
								println( "Average distance: " + avgDistance );
								println( "Hausdorff distance: " + hausdorffDistance + "\n" );
	
						} );
				
						writeModelFile( posterioriModel , tempLoc + "/partialModel_" + targetNo + ".h5" )
						}

						val posteriorGroup = UI.createGroup( "Posterior "+ targetNo );
						UI.show( posteriorGroup, posterioriModel, "StatisticalModel" );
			    	
						//return most likely sample as result
						posterioriModel.mean
					
		 		} ).toArray;

			 return meshesCompleted;					      

		}
		
}

}