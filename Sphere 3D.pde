// homogenous 3D vector
class hom3d {
  float x1, x2, x3, x4;

  hom3d(float x1, float x2, float x3, float x4) {
    this.x1=x1;
    this.x2=x2;
    this.x3=x3;
    this.x4=x4;
  }
}

// inhomogenous 3D vector
class inhom3d {
  float x, y, z;

  inhom3d(float x, float y, float z ) {
    this.x=x;
    this.y=y;
    this.z=z;
  }
}

// =============================================================================

// inhomogenous => homogenous conversion
hom3d toHom3d( inhom3d A ) {
  return new hom3d( A.x, A.y, A.z, 1 );
}

// homogenous => inhomogenous conversion
inhom3d toInhom3d( hom3d a ) {
  return new inhom3d( a.x1 / a.x4, a.x2 / a.x4, a.x3 / a.x4 );
}

// matrix * vector; returns M * P
hom3d mulPoint( float M[][], hom3d P ) {
  hom3d P_ = new hom3d(  M[0][0] * P.x1 + M[0][1] * P.x2 + M[0][2] * P.x3 + M[0][3] * P.x4, 
    M[1][0] * P.x1 + M[1][1] * P.x2 + M[1][2] * P.x3 + M[1][3] * P.x4, 
    M[2][0] * P.x1 + M[2][1] * P.x2 + M[2][2] * P.x3 + M[2][3] * P.x4, 
    M[3][0] * P.x1 + M[3][1] * P.x2 + M[3][2] * P.x3 + M[3][3] * P.x4);
  return P_;
}

// matrix * matrix; C = A * B (stores A * B in C)
void mulMatrix( float A[][], float B[][], float C[][] ) {
  float sum;
  for (int i = 0; i <= 3; i++) {
    for (int j = 0; j <= 3; j++) {
      sum = 0;
      for (int k = 0; k <= 3; k++) {
        sum = sum + A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}

// returns the vector pointing from A to B (B - A)
inhom3d vectorFromTo( inhom3d A, inhom3d B ) {
  return new inhom3d ( B.x - A.x, B.y - A.y, B.z - A.z  );
}

// returns the cross product of a and b (a x b)
inhom3d crossProd( inhom3d a, inhom3d b ) {
  return new inhom3d ( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x );
}

// returns the dot product of a and b (a * b)
float dotProd( inhom3d a, inhom3d b ) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

// returns the length of the vector (||a|| = sqrt(a * a))
float vectorLength( inhom3d a ) {
  return sqrt(dotProd(a, a));
}

// returns the vector normalized to unit length (a / ||a|| or a / length(a))
inhom3d normalizeVector( inhom3d a ) {
  float l = vectorLength(a);
  return new inhom3d( a.x / l, a.y / l, a.z / l );
}

// evaluates the coordinate function of the sphere at the parameters phi and theta
inhom3d spherePoint(float phi, float theta) {
  return new inhom3d( 
    cos(theta) * sin(phi),
    sin(theta) * sin(phi),
    cos(phi));
}

// =============================================================================
// number of sides the sphere has (in both directions)
int sphereSides = 32;

// phi in [0, PI] => step/delta = PI / sides
float phiStep = PI / sphereSides;

// theta in [0, 2PI] => step/delta = 2PI / sides
float thetaStep = (2 * PI) / sphereSides; 

// center of projection
float s = 150;
inhom3d center = new inhom3d( 0, 0, s );

// normalized vector pointing towards the light source
inhom3d toLight = normalizeVector( new inhom3d( 100, 0, 20 ) );

// =============================================================================
// rotation values
float alphax = 60;
float alphay = 50;

// translation matrix
float [][] matrixT = { 
  {1, 0, 0, 512}, 
  {0, 1, 0, 384}, 
  {0, 0, 1, 0}, 
  {0, 0, 0, 1}  };

// scale matrix
float [][] matrixS = { 
  {125, 0, 0, 0}, 
  {0, 125, 0, 0}, 
  {0, 0, 125, 0}, 
  {0, 0, 0, 1}  };

// projection matrix
float [][] matrixP = { 
  {1, 0, 0, 0}, 
  {0, 1, 0, 0}, 
  {0, 0, 0, 0}, 
  {0, 0, -1.0/s, 1}  };

// matrices to hold the temporary values while computing the final sum
float [][] matrixTmp1 = new float [4][4];
float [][] matrixTmp2 = new float [4][4];
float [][] matrixTmp3 = new float [4][4];
float [][] matrixTmp4 = new float [4][4];
float [][] matrixTmp5 = new float [4][4];

// final composite transformation matrices
float [][] matrixM3D = new float [4][4];
float [][] matrixM = new float [4][4];

// =============================================================================

void setup()
{
  size( 1024, 768 );
  ellipseMode( RADIUS );
}

// =============================================================================

void draw() {

  background(200);

  // rotation matrix (x-axis)
  float [][] matrixRx = { 
    { 1, 0, 0, 0}, 
    {0, cos(alphax), -sin(alphax), 0}, 
    {0, sin(alphax), cos(alphax), 0}, 
    {0, 0, 0, 1}  };

  // rotation matrix (y-axis)
  float [][] matrixRy = { 
    { cos(alphay), 0, sin(alphay), 0}, 
    {0, 1, 0, 0}, 
    {-sin(alphay), 0, cos(alphay), 0}, 
    {0, 0, 0, 1}  };

  // compute the necessary transformation matrices here
  // TODO
  
  //M3D = Rx * Ry
  mulMatrix(matrixRx, matrixRy, matrixTmp1);
  //M3D = Rx * Ry * S
  mulMatrix(matrixS, matrixTmp1, matrixM3D);
  //M = P * M3D
  mulMatrix(matrixP, matrixM3D, matrixTmp2);
  //M = T * P * M3D
  mulMatrix(matrixT, matrixTmp2, matrixM);

  // go through all sides of the sphere
  for (int i = 0; i < sphereSides; i++ ) {
  for (int j = 0; j < sphereSides; ++j) {
  	// phi and theta parameters for the current face
  	float phi = i * phiStep;
  	float theta = j * thetaStep;
  	float phiNext = (i + 1) * phiStep;
  	float thetaNext = (j + 1) * thetaStep;
  	
  	// compute the 4 vertices of the face
    // TODO
    inhom3d p0, p1, p2, p3;
    
    p0 = (spherePoint(phi,theta));
    p1 = (spherePoint(phiNext,theta));
    p2 = (spherePoint(phiNext,thetaNext));
    p3 = (spherePoint(phi,thetaNext));

    // coordinates of the rotated and scaled vertices of the face
    // TODO
    
    //P'i = M3D * Pi
    inhom3d rotatedP0 = toInhom3d(mulPoint(matrixM3D, toHom3d(p0)));
    inhom3d rotatedP1 = toInhom3d(mulPoint(matrixM3D, toHom3d(p1)));
    inhom3d rotatedP2 = toInhom3d(mulPoint(matrixM3D, toHom3d(p2)));
    inhom3d rotatedP3 = toInhom3d(mulPoint(matrixM3D, toHom3d(p3)));

    //edges of the transformed face
    inhom3d e0 = vectorFromTo(rotatedP0, rotatedP1);
    inhom3d e1 = vectorFromTo(rotatedP0, rotatedP2);
    inhom3d e2 = vectorFromTo(rotatedP0, rotatedP3);

    // normal vector of the transformed face
	  // TODO
    inhom3d normal;
    if (i < sphereSides - 1)
	     //normal = (e0 x e1) / || e0 x e1 ||
       normal = normalizeVector(crossProd(e0,e1))
	    ;
	  else
	     //normal = (e0 x e2) / || e0 x e2 ||
       normal = normalizeVector(crossProd(e0,e2))
	    ;
    
    // amount of light illuminating the face
    // TODO
    //Shading
    float intensity = (dotProd(normal, toLight) + 1.0f) / 2.0f;
    color col = color(intensity * 200, intensity * 0, intensity * 0);

    // skip the face if not visible
    // TODO
    //Back-Face culling
	  boolean visible = dotProd(normal, vectorFromTo(rotatedP0, center)) >= 0;
    if (!visible)
      continue;
    
    // projected vertices
    //P''i = M * Pi
    inhom3d projectedP0 = toInhom3d(mulPoint(matrixM, toHom3d(p0)));
    inhom3d projectedP1 = toInhom3d(mulPoint(matrixM, toHom3d(p1)));
    inhom3d projectedP2 = toInhom3d(mulPoint(matrixM, toHom3d(p2)));
    inhom3d projectedP3 = toInhom3d(mulPoint(matrixM, toHom3d(p3)));

    //draw the face
    fill(col);
    beginShape();
    vertex(projectedP0.x, projectedP0.y);
    vertex(projectedP1.x, projectedP1.y);
    vertex(projectedP2.x, projectedP2.y);
    vertex(projectedP3.x, projectedP3.y);
    endShape(CLOSE);
  }
}
  // keep rotating
  alphax += PI / 360;
  alphay += PI / 360;
}
