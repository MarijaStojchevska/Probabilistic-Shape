package example;

import java.io._

import scalismo.ui.api.ScalismoUI
import scalismo.ui.api.{ Group ⇒ apiGroup }

import scalismo.common._
import scalismo.registration._
import scalismo.io._
import scalismo.geometry._
import Array._
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.dataset.DataCollection

//import scala.collection.parallel.ParIterableLike.Reduce

package AlignmentAndCompletionTrials {

object AlignmentTrials {
	
		/**
		 * wrapper to store a mesh and landmark tuple
		 */
		case class MeshLandmark( mesh : TriangleMesh[ _3D ] , 
														 lm		: Seq[ Landmark[ _3D ] ] ){}
		
		/**
		 * wrapper to store a mesh array and a landmark array as tuple
		 */
		case class MeshLandmarkArrays( meshes : Array[ TriangleMesh[ _3D ] ] , 
																	 lms 		: Array[ Seq[ Landmark[ _3D ] ] ] ){}
	
		/**
		 * read all meshes from location 
		 * ! here, we assume all local files to be meshes
		 * returns null, if files do not exist !
		 */
		def readMeshes( location : String ) : Array[ TriangleMesh[ _3D ] ] = {
			
				val files = new File( location ).listFiles();
				
				if( files.isEmpty ) return null
				
				val meshes = files.map( file => {
					
						val mesh = MeshIO.readMesh( file ).get;
						
						mesh
					
				});
				
				return meshes;
			
		}

		/**
		 * read specific mesh by given path
		 * returns null, if file does not exist !
		 */
		def readMesh( name : String ) : TriangleMesh[ _3D ] = {
			
				val file = new java.io.File( name );
				if( !file.exists() ) return null;

				val mesh = MeshIO.readMesh( file ).get;				
				
				return mesh;
			
		}

		/**
		 * read specific landmark by given path
		 * returns null, if file does not exist !
		 */
		def readLM( name : String ) : Seq[ Landmark[ _3D ] ] = {

				val file = new java.io.File( name );
				if( !file.exists() ) return null;

				val lm = LandmarkIO.readLandmarksJson[ _3D ]( file ).get;
				
				return lm;
			
		}

		/**
		 * read specific mesh & landmark by given path
		 * returns null, if file does not exist !
		 */
		def readMeshLM( targetMesh : String , targetLM : String ) : MeshLandmark = {
			
				val fileMesh = new File( targetMesh );
				if( !fileMesh.exists() ) return MeshLandmark( null , null );
				
				val mesh = MeshIO.readMesh( fileMesh ).get;
				
				val fileLM = new File( targetLM );
				if( !fileLM.exists() ) return MeshLandmark( null, null );
				
				val lm = LandmarkIO.readLandmarksJson[ _3D ]( fileLM ).get;
			
				return MeshLandmark( mesh , lm );
				
		}
	
		/**
		 * read all meshes and landmarks from given locations
		 * ! here, we assume all local files to be meshes
		 * returns null, if file do not exist !
		 */
		def readMeshsLMs( locationMeshes : String , locationLMs : String ) : MeshLandmarkArrays = {
			
				val fileMeshes = new File( locationMeshes ).listFiles();
				if( fileMeshes.isEmpty ) return ( MeshLandmarkArrays( null , null ) );
				
				val meshes : Array[ TriangleMesh[ _3D ] ] = 
																					fileMeshes.map( file ⇒  MeshIO.readMesh( file ).get );

				val fileLMs = new File( locationLMs ).listFiles();
				if( fileLMs.isEmpty ) return ( MeshLandmarkArrays( null , null ) );

				val LMs : Array[ Seq[ Landmark[ _3D ] ] ] = fileLMs.map( 
																						file ⇒  LandmarkIO.readLandmarksJson[ _3D ]( file ).get );

				return MeshLandmarkArrays( meshes , LMs );
				
		}
		
		/**
		 * write tuple of mesh- und landmark-arrays to given locations
		 * ! here, we assume the locations exists and existing files can be overwritten
		 */
		def writeMeshsLMs( locationMeshes : String , locationLMs : String , data : MeshLandmarkArrays ) {
			
				data.meshes.map( mesh ⇒  MeshIO.writeMesh( mesh , new java.io.File( locationMeshes ) ).get );

				data.lms.map( lm ⇒  LandmarkIO.writeLandmarksJson( lm , new java.io.File( locationMeshes ) ) );
				
		}

		/**
		 * use fully automated pre-implemented Procrast Data Collection read
		 * the resulting mesh data is generated from single reference mesh and the deformation fields
		 * describing all other meshes
		 * the landmarks are read manually and aligned using the before mentioned deformation fields
		 * returns null, if file does not exist !
		 */
		def readMeshsLMsAndAlign( locationMeshes : String  , 
														 	locationLMs : String	) : ( MeshLandmark , MeshLandmarkArrays ) = {

				implicit val rng = scalismo.utils.Random(42);
			
				val fileMeshes = new File( locationMeshes ).listFiles();
				if( fileMeshes.isEmpty ) return null;

				val directory  = new File( locationMeshes );
				
				// reading first mesh as reference
				val meshFirstRef : TriangleMesh[ _3D ] = MeshIO.readMesh( fileMeshes( 0 ) ).get;
								
				val fileLMs = new File( locationLMs ).listFiles();
				
				var LMs : Array[ Seq[ Landmark[ _3D ] ] ] = fileLMs.map( 
																						file ⇒  LandmarkIO.readLandmarksJson[ _3D ]( file ).get );

				// primingn Procrustes Analysis to given mesh dataset ( by folder )
				val GenProcrastDataObj = DataCollection.fromMeshDirectory( meshFirstRef , directory )._1.get;
				
				// choosing upper iteration length of 20 and a goal precission of 0.1
				val GenProcrastDataCollection = DataCollection.gpa( GenProcrastDataObj , 20 , 0.1 )

				// return the mean transformation to generate a mean reference landmark set
				val GenProcrastTransf = GenProcrastDataCollection.meanTransformation;
				
				var theMean : MeshLandmark = MeshLandmark( null , null );
				// return the mean shape as reference
				var meanMesh = GenProcrastDataCollection.meanSurface;

				// generate the training set using the transformations for the mean
				// ! this will generate an interpolated version of the training data !
				val meshes = GenProcrastDataCollection.dataItems.seq.map( dataelem ⇒ {
					
					meanMesh.transform( dataelem.transformation );
					
				} ).toArray

				// generate the mean landmark set using the transformation from the original 
				// Procrustes reference
				val meanLM = ( 0 until LMs( 0 ).size ).map( elem ⇒ {
					
					LMs( 0 )( elem ).transform( GenProcrastTransf );
					
				});
				
				// generate the landmarks analogue to the meshes from the mean LM
				LMs = ( 0 until GenProcrastDataCollection.dataItems.seq.size ).map( dataelem ⇒ {
					
					( 0 until LMs( dataelem ).size ).map( lm ⇒ {
						
							meanLM( lm ).transform( 
																GenProcrastDataCollection.dataItems.seq( dataelem ).transformation );
						
					});
					
				} ).toArray


				return ( MeshLandmark( meanMesh , meanLM ) , MeshLandmarkArrays( meshes , LMs ) );
				
		}

		/**
		 * calculates mean mesh with a fitting set of landmarks using mean deformation
		 * this is a non-iterative version :
		 * 		it simply calculate the mean of a given mesh-landmark dataset ( non ICP )
		 */
		def getMean( data : ⇒ MeshLandmarkArrays ) : MeshLandmark = {
			
	      // setting a seed for the random generator to allow for reproducible results
	      implicit val rng = scalismo.utils.Random(42)
	      
	      // using the first mesh as a reference, average over all deformation (differences to the ref)
				val meanDeformation = data.meshes( 0 ).pointSet.pointsWithId.map( pt ⇒ {
						
						var meanPtDeformation = EuclideanVector( 0 , 0 , 0 );
						data.meshes.map( mesh ⇒ {
								
								val refpt = mesh.pointSet.findClosestPoint( pt._1 );
								val diff  = refpt.point - pt._1;
								meanPtDeformation += diff / data.meshes.size;
							
						});
						meanPtDeformation
						
				} ).toIndexedSeq;

	      // create continuous deformation
	      
				val deformFieldDiscrete = DiscreteField[ _3D , UnstructuredPointsDomain[ _3D ], 
																								 EuclideanVector[ _3D ] ]( data.meshes( 0 ).pointSet , 
																										 											 meanDeformation );

		    val interpolator 				= NearestNeighborInterpolator[ _3D , EuclideanVector[ _3D ] ]();
		    val deformFieldContin   = deformFieldDiscrete.interpolate( interpolator );
		    val transformation			= Transformation( ( pt : Point[ _3D ] ) ⇒  
		    																															pt + deformFieldContin( pt ) );

				// create mean set using the initial reference
		    val meanMesh = data.meshes( 0 ).transform( transformation );
				val meanLM   = ( 0 until data.lms( 0 ).size ).map( elem ⇒ {
																							data.lms( 0 )( elem ).transform( transformation ) } );

				return MeshLandmark( meanMesh , meanLM );
			
		}

		/**
		 * calculates a mean set and uses it as reference aligning all other data sets to it
		 * ( single run version ... we presume a certain realignment to generate a meaningful mean )
		 * Given landmarks are used as reference ( non ICP )
		 */
		def getMeanAndAlign( data : ⇒ MeshLandmarkArrays ) : MeshLandmark = {

				// calculate the mean
				val mean = getMean( data );

				// align all dataset to gotten mean
				alignMeshesAndLandmarksToReference( data , mean );
				
				return mean;
			
		}

		/**
		 * uses the reference landmark set to generate landmark sets for all given meshes
		 * ( this is no IPC version , but a simple Closest Point assumption )
		 */
		def alignLandmarksTo( referenceLMs : Seq[ Landmark[ _3D ] ] , 
													meshes : ⇒ Array[ TriangleMesh[ _3D ] ] ) : 
																																Array[ Seq[ Landmark[ _3D ] ] ] = {
			
				var setOfNewLMs = new Array[ Seq[ Landmark[ _3D ] ] ]( meshes.size );
				
				setOfNewLMs = meshes.map( mesh ⇒ {
					
						val lmSeq : Seq[ Landmark[ _3D ] ] = ( 0 until referenceLMs.size ).map( elem ⇒ {
	
								val point				 = referenceLMs( elem ).point;
								val pointClosest = mesh.pointSet.findClosestPoint( point );
								val lm : Landmark[ _3D ] = Landmark( referenceLMs( elem ).id , pointClosest.point );
								lm
								
						} );						
						lmSeq
				
				} );
				return setOfNewLMs;
				
		}

		/**
		 * Iterative alignment by calculating a mean set and align to it. 
		 * Given landmarks are used as reference ( non ICP )
		 * Repeat until : 
		 * 			- new found mean shows a total square norm difference to the previous one below 0.1
		 * 			- maximum number of given iterations is reached
		 */
		def alignMeshesAndLandmarks( data : ⇒  MeshLandmarkArrays , iterations : Int ) : MeshLandmark = {

				implicit val rng = scalismo.utils.Random(42);

				// variable initialisation
				var diff : Double 	= 1.0;
				var counter : Int   = iterations;
				var mean						= MeshLandmark( data.meshes( 0 ) , data.lms( 0 ) ); // the reference set
				
				//Procrustes Alignment manually
				while( diff > 0.1 && counter > 0 ){
					
						//regid alignemnt using provided landmarks
						val rigidTransformation = data.lms.map( lm ⇒ {
							
								LandmarkRegistration.rigid3DLandmarkRegistration( lm , mean.lm , 
																																	Point3D( 0 , 0 , 0 ) )
							
						} );

						// align landmarks
						( 0 until data.lms.size ).map{ elem ⇒
							
								data.lms( elem ) = data.lms( elem ).map( lm ⇒ 
																											lm.transform( rigidTransformation( elem ) ) );
						
						}

						// align meshes
						( 0 until data.meshes.size ).map{ elem ⇒
							
								data.meshes( elem ) = data.meshes( elem ).transform( rigidTransformation( elem ) );
							
						}

						// calculate new mean and align to it
						val newMean = getMeanAndAlign( data )

						// evaluate new mean based on previous one
						diff = 0;
						( 0 until newMean.mesh.pointSet.numberOfPoints ).map( elem ⇒ {
							
								diff = diff + ( newMean.mesh.pointSet.points.toSeq( elem ) - 
																mean.mesh.pointSet.points.toSeq( elem ) ).norm2;
							
						} );
						
						// choose quality of mean-ing
						counter  = counter - 1;
						if( diff > 0.1 && counter > 0 ) {
							
							mean = newMean;
						
						}

				}
				
				return mean;
			
		}
		
		/**
		 * Non-iterative alignment by aligning the data to calculating a mean set and align to it. 
		 * Given landmarks are used as reference ( non ICP )
		 */
		def alignMeshesAndLandmarks_singleRun( data : ⇒ MeshLandmarkArrays	) : MeshLandmark = {
		
				implicit val rng = scalismo.utils.Random(42);
			
				var meanMesh				= data.meshes( 0 );
				var meanLMs					= data.lms( 0 );
				
				//regid alignemnt manually...
				
				// calcualte transformations
				val rigidTransformation = data.lms.map( lm ⇒ {
							
						LandmarkRegistration.rigid3DLandmarkRegistration( lm , meanLMs , 
																															Point3D( 0 , 0 , 0 ) )
					
				} );
				
				// align landmarks
				( 0 until data.lms.size ).map{ elem ⇒
					
						data.lms( elem ) = data.lms( elem ).map( lm ⇒ 
																									lm.transform( rigidTransformation( elem ) ) );

				}

				// align meshes
				( 0 until data.meshes.size ).map{ elem ⇒
					
						data.meshes( elem ) = data.meshes( elem ).transform( rigidTransformation( elem ) );

				}

				// calculate mean and align to it
				val mean = getMeanAndAlign( data );

				return mean;
			
		}
		
		/**
		 * align data set by given reference set using landmarks as points of reference ( non ICP )
		 */
		def alignMeshesAndLandmarksToReference( newData : ⇒ MeshLandmarkArrays , 
																						reference : ⇒ MeshLandmark ){
			
				implicit val rng = scalismo.utils.Random(42);
				
				//regid alignemnt manually...
				
				// calculate transformations
				val rigidTransformation = newData.lms.map( lm ⇒ {
							
						LandmarkRegistration.rigid3DLandmarkRegistration( lm , reference.lm , 
																															Point3D( 0 , 0 , 0 ) )
					
				} );

				// align landmarks
				( 0 until newData.lms.size ).map{ elem ⇒
					
						newData.lms( elem ) = newData.lms( elem ).map( lm ⇒ 
																									lm.transform( rigidTransformation( elem ) ) );
				
				}

				// align meshes
				( 0 until newData.meshes.size ).map{ elem ⇒
					
						newData.meshes( elem ) = newData.meshes( elem ).transform( rigidTransformation( elem ) );
					
				}
			
		}
		
		/**
		 * Calculates the centreoids of all data sets using all mesh points. Using these centreoids as
		 * reference, the data sets are aligned to the first data set
		 */
		def centeroidAlignment( data : ⇒ MeshLandmarkArrays	) {
																 						 
				implicit val rng = scalismo.utils.Random(42);

				// calculate centreoids
				val centeroids = data.meshes.map( mesh ⇒ {
					
			      val vector = mesh.pointSet.points.toSeq.map{point => point.toVector}
			      val total_mesh_points = mesh.pointSet.numberOfPoints
			      (vector.reduce((vec1, vec2) => vec1 + vec2)) / total_mesh_points
					
				} );

				//build translation tranformation based on centreoids ( use first one as reference )
				val transf = centeroids.map( ptSet ⇒  TranslationTransform[ _3D ]( ptSet.*( -1 ) ) );

				// align meshes
				( 0 until data.meshes.size ).map( elem ⇒
						
						data.meshes( elem ) = data.meshes( elem ).transform( transf( elem ) )
				
				) 

				// align landmarks
				( 0 until data.lms.size ).map{ elem ⇒
					
						data.lms( elem ) = data.lms( elem ).map( lm ⇒  lm.transform( transf( elem ) ) );
				
				}

		} 

		/**
		 * depricated
		 */
		def readFilesAndPreAlign() : MeshLandmarkArrays = {
			
				implicit val rng = scalismo.utils.Random(42);
			
				val numberOfElements : Int = 50;
				
				// load all meshes and according landmarks and align them to one another 
				// (using 1st as reference)
	      var meshes     = new Array[ TriangleMesh[ _3D ] ]( numberOfElements );
	      var landmarks  = new Array[ Seq[ Landmark[ _3D ] ] ]( numberOfElements );
	      meshes( 0 ) 	 = MeshIO.readMesh(new File( s"data/meshes/${ 0 }.stl" ) ).get;
	      landmarks( 0 ) = LandmarkIO.readLandmarksJson[ _3D ]( 
	      																						new File( s"data/landmarks/${ 0 }.json" ) ).get;
	      for( elem ← 1 until numberOfElements by 1 ){
	
	      		// reading from files
		      	meshes( elem ) 		= MeshIO.readMesh( new File( s"data/meshes/${ elem }.stl" ) ).get;
		      	landmarks( elem ) = LandmarkIO.readLandmarksJson[ _3D ]( new File( 
	      																									s"data/landmarks/${ elem }.json" ) ).get;
		      	
		      	// getting transformation for alignment
		      	var transf : RigidTransformation[_3D] = LandmarkRegistration.
		      													rigid3DLandmarkRegistration( landmarks( elem ) , landmarks( 0 ) , 
		      																											 Point( 0 , 0 , 0 ) );

		      	// applying it to mesh and landmarks to align input samples
		      	meshes( elem )    = meshes( elem ).transform( transf );
		      	var transfLM = new Array[ Landmark[ _3D ] ]( landmarks( 0 ).size );
				    for( subElem ← 0 until landmarks( 0 ).size ){
				    	
				    		transfLM( subElem ) = Landmark( landmarks( 0 )( subElem ).id , 
		    																		landmarks( elem )( subElem ).transform( transf ).point );
				    	
				    }
		      	landmarks( elem ) = transfLM;
	
	      }
	      
	      return MeshLandmarkArrays( meshes , landmarks );
			
		}
		
		/**
		 * depricated
		 */
		def meanAlignmentVariance( meshes : ⇒  Array[ TriangleMesh[ _3D ] ] , 
															 landmarks : ⇒ Array[ Seq[ Landmark[ _3D ] ] ] , 
															 ui : ⇒  ScalismoUI , uiGroup : ⇒ apiGroup ) : 
															 							( Seq[ Landmark[ _3D ] ] , Array[ Array[ Float ] ] ) = {
			
				implicit val rng = scalismo.utils.Random(42);

	      val numberOfElements : Int = meshes.length;

	      var meanLM = new Array[ Landmark[ _3D ] ]( landmarks( 0 ).size );

			      //************************************************************************************//
			      // calculate mean shape
		
			      var meanVect = new Array[ EuclideanVector[ _3D ] ]( landmarks( 0 ).size );
				    for( subElem ← 0 until landmarks( 0 ).size by 1 ){
				      	
				   	 		meanVect( subElem ) = EuclideanVector( 0.0 , 0.0 , 0.0 );
				      	
				    } 
			      for( elem ← 0 until numberOfElements by 1 ){
			      	
			      	
					      for( subElem ← 0 until landmarks( 0 ).size by 1 ){
					      	
					      		meanVect( subElem ) = meanVect( subElem ) + 
					      												  landmarks( elem )( subElem ).point.toVector;
					      	
					      }
			      	
			      }
				    for( subElem ← 0 until landmarks( 0 ).size by 1 ){
				      	
				    		meanVect( subElem ) = meanVect( subElem ) / landmarks.size;
				      	
				    }
				    
				    // array to Landmark/Seq
				    for( subElem ← 0 until meanLM.size by 1 ){
				    	
				    		meanLM( subElem ) = Landmark( landmarks( 0 )( subElem ).id , 
				    																	meanVect( subElem ).toPoint )
				    	
				    }

			      //************************************************************************************//      
			      // realign shapes to mean
		
				    for( elem ← 0 until numberOfElements by 1 ){
			
					      	var transf 		 	  = LandmarkRegistration.rigid3DLandmarkRegistration( 
					      																	landmarks( elem ) , meanLM , Point( 0 , 0 , 0 ) );
					      	meshes( elem ) 	  = meshes( elem ).transform( transf );
					      	var transfLM = new Array[ Landmark[ _3D ] ]( landmarks( 0 ).size );
							    for( subElem ← 0 until transfLM.size by 1 ){
							    	
								    	transfLM( subElem ) = Landmark( landmarks( 0 )( subElem ).id , 
								    												landmarks( elem )( subElem ).transform( transf ).point );
							    	
							    }
			      	
			      }
				    
				    				    
//			      //************************************************************************************//      
//			      // evaluate effect
//		
//				    val meanLM_old  = listOFMeans.last;
//				    listOFMeans 	  = listOFMeans + meanLM;
//				    
//				    sum = 0.0;
//				    for( lm ← 0 until meanLM.size by 1 ){
//				    	
//				    		val delta : EuclideanVector[ _3D ] = meanLM( lm ).point - meanLM_old( lm ).point;
//				    		val dist : Double = delta.norm;
//				    		sum = sum + dist;
//				    	
//				    }
//			    
//	      }
//	      ui.show( uiGroup , listOFMeans.toSeq.flatten , "femurMean" );
      
	      //****************************************************************************************//
	      // calculating co-variance matrix
	      
		    // covariant matrix is of size 9 * n^2 with n = # landmarks to be considered
				val dim_in_blocks  : Int = landmarks( 0 ).size;
				val size_in_blocks : Int = dim_in_blocks * dim_in_blocks;
				val dim_block      : Int = 3;
				val size_block 		 : Int = dim_block * dim_block;
				val size           : Int = size_in_blocks * size_block;
				val dim						 : Int = dim_in_blocks * dim_block;
				val numSam				 : Int = landmarks.length;
				
				var covari = ofDim[ Float ]( dim , dim ); // use correct size
				// the covariant matrix looks like : for all i,j dimensions and k,l landmarks :  
				//    ( x_i_k - mean_k ) * ( x_j_k - mean_k ) , ( x_i_k - mean_k ) * ( x_j_l - mean_l )
				//    ( x_j_k - mean_k ) * ( x_i_l - mean_l ) , ( x_j_k - mean_k ) * ( x_j_k - mean_k )
				// each representing a lock with any i,j and k,l combination
				for( run_i ← 0 until dim by 1 ){
					
						for( run_j ← 0 until dim by 1 ){
							
							val lm_i  = run_i / dim_block;
							val dim_i = run_i % dim_block;
							val lm_j  = run_j / dim_block;
							val dim_j = run_j % dim_block;
							var sum : Float = 0;
							for( sample ← 0 until numSam by 1 ){
								
									sum = sum + ( ( landmarks( sample )( lm_i ).point( dim_i ).toDouble - 
																	meanLM( lm_j ).point( dim_i ).toDouble ) * 
																	( landmarks( sample )( lm_i ).point( dim_j ).toDouble - 
																	meanLM( lm_j ).point( dim_j ).toDouble ) ).toFloat;
								
							}
							sum 										 = sum / numSam;
							covari( run_i )( run_j ) = sum;
											
						}
					
				}
				
				return ( meanLM , covari );
			
		}
		
		/**
		 * depricated
		 */
		def writeMeanCoVariance( meanLM : Seq[ Landmark[ _3D ] ] , covariLM : Array[ Array[ Float ] ] ){
			
	      val writer = new PrintWriter( new File( "MeanAndVariance.txt" ) );
	      writer.write( "mean landmarks : \n" );
	      for( lm ← meanLM ) {
	      	
	      		writer.write( "id=" );
	      		writer.write( lm.id );
	      		writer.write( " , " );
	      	  writer.write( lm.point.toString() );
	      	  writer.write( " } ; " )
	      	
	      }
	      writer.write( "\n Co-Variance : \n" );
	      writer.write( "{ " );
	      for( row ← 0 until covariLM.size by 1 ){
	      	
	      		if( row > 0 ) writer.write( " ; " );
	      		for( column ← 0 until covariLM( 0 ).size by 1 ){
	      			
	      				if( column > 0 ) writer.write( " , " );
	      				writer.write( "" + covariLM( row )( column ) );
	      			
	      		}

	      }
	      writer.write( " }" );
      	writer.close();
			
		}
		
		/**
		 * depricated
		 */
		def aligningAndAnalysis( ui : ⇒  ScalismoUI ) {

	      // setting a seed for the random generator to allow for reproducible results
	      implicit val rng = scalismo.utils.Random(42)
	
	      // read a mesh from file
	      val femurGroup = ui.createGroup( "Femur Complition" );

	      // reading & Pre-Aligning Meshes & Landmarks
	      var meshesAndLandmarks = readFilesAndPreAlign();

	      val numberOfElements : Int = meshesAndLandmarks.meshes.length;
	      
	      val ( meanLM , covariLM ) = meanAlignmentVariance( meshesAndLandmarks.meshes , 
	      																									 meshesAndLandmarks.lms , ui , femurGroup );

				// printing the gained values
				
//	      var alignedMeshViews = new Array( numberOfElements );
//	      var alignedLMViews = new Array( numberOfElements );
//	      for( elem ← 0 until numberOfElements by 1 ){
//	      	
//	      	alignedMeshViews( elem ) = ui.show( femurGroup , meshes( elem ) , s"femur${elem}}" );
//	      	alignedLMViews( elem )   = ui.show( femurGroup , landmarks( elem ) , s"femur${elem}}" );
//	      	
//	      }
	      var muLM : Seq[ Landmark[ _3D ] ] = meanLM;
	      ui.show( femurGroup , muLM , "femurMean" );

	      writeMeanCoVariance( muLM , covariLM );
      	
  	}

}

}
