#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

// input : 
// h , w , msgLen-> submatrix and matrix
// x , msg -> y
//
// eg:( no blank space in the end of each line ):
// 2 /h
// 2 /w
// 4 /msgLen
// 1 0 
// 1 1 /submatrix
// 1 0 1 1 0 0 0 1 /x
// 0 1 1 1 /m


int h = 0;
int w = 0;
int msgLen = 0;

int numState = 0;
int pathLen = 0;

unsigned int *x = NULL;
unsigned int *y = NULL;
unsigned int *msg = NULL;

unsigned int **subMatrix = NULL;
unsigned int **matrix = NULL;
unsigned int **Hmatrix = NULL;

double *wght = NULL;
double *newwght = NULL;

unsigned int indx = 0;
unsigned int indm = 0;

unsigned int **path = NULL;

void readInput( char *filename )
{
	FILE *file = fopen( filename , "r" );
	fscanf(file , "%d\n%d\n%d\n" , &h , &w , &msgLen);
	numState = pow(2 , h);
	pathLen = w * msgLen;

	int i = 0;

	subMatrix = (unsigned int **)calloc( h , sizeof( unsigned int * ) );
	for( i = 0 ; i < h ; i++ )
		subMatrix[i] = (unsigned int *)calloc( w , sizeof( unsigned int ) );

	int j = 0;
	for( i = 0 ; i < h ; i++ )
	{
		for( j = 0 ; j < w ; j++ )
		{
			subMatrix[i][j] = fgetc(file) - '0';
			fgetc(file);
		}
	}

	printf("subM :\n");
	for( i = 0 ; i < h ; i++ )
	{
		for( j = 0 ; j < w ; j++ )
			printf("%d " , subMatrix[i][j]);
		printf("\n");
	}

	msg = (unsigned int *)calloc( msgLen , sizeof( unsigned int ) );
	y = (unsigned int *)calloc( pathLen , sizeof( unsigned int ) );
	x = (unsigned int *)calloc( pathLen , sizeof( unsigned int ) );

	for( i = 0 ; i < pathLen ; i++ )
	{
		x[i] = fgetc(file) - '0';
		fgetc(file);
	}

	for( i = 0 ; i < msgLen ; i++ )
	{
		msg[i] = fgetc(file) - '0';
		fgetc(file);
	}

	fclose( file );
}

void init()
{
	int i = 0;
	int j = 0;

	unsigned int columnVal[w];
	for( i = 0 ; i < w ; i++ )
	{
		unsigned int tmp = 0;
		for( j = h - 1 ; j >= 0 ; j-- )
		{
			tmp <<= 1;
			tmp += subMatrix[j][i];
		}
		columnVal[i] = tmp;
		printf("columnVal[%d] = %d\n" , i , columnVal[i]);
	}

	matrix = (unsigned int**)malloc( msgLen * sizeof( unsigned int * ) );
	for( i = 0 ; i < msgLen ; i++ )
		matrix[i] = (unsigned int *)malloc( w * sizeof( unsigned int ) );
	for( i = 0 ; i < msgLen ; i++ )
	{
		unsigned int rem = pow(2 , msgLen - i );
		for( j = 0 ; j < w ; j++ )
			matrix[i][j] = columnVal[j] % rem;
	}

	printf("matrix :\n");
	for( i = 0 ; i < msgLen ; i++ )
	{
		for( j = 0 ; j < w ; j++ )
			printf("%d " , matrix[i][j]);
		printf("\n");
	}

	path = (unsigned int **)malloc( numState * sizeof( unsigned int * ) );
	for( i = 0 ; i < numState ; i++ )
		path[i] = (unsigned int *)calloc( pathLen , sizeof( unsigned int ) );

	Hmatrix = (unsigned int **)malloc( msgLen * sizeof( unsigned int * ) );
	for( i = 0 ; i < msgLen ; i++ )
		Hmatrix[i] = (unsigned int *)calloc( pathLen , sizeof( unsigned int ) );

	for( j = 0 ; j < pathLen ; j++ )
	{
		unsigned int delta = j / w;
		for( i = 0 ; i < h ; i++ )
		{
			if( i + delta >= msgLen )
				break;
			Hmatrix[ i + delta ][ j ] = subMatrix[ i ][ j % w ]; 	
		}
	}

	printf("matrix H:\n");
	for( i = 0 ; i < msgLen ; i++ )
	{
		for( j = 0 ; j < pathLen ; j++ )
			printf("%d " , Hmatrix[i][j]);
		printf("\n");
	}

	newwght = (double *)calloc( numState , sizeof( double ) );
	wght = (double *)calloc( numState , sizeof( double ) );
	wght[0] = 0.0;
	for( i = 1 ; i < numState ; i++ )
		wght[i] = INFINITY;

}

double distortion()
{
	// distortion function
	return 1.0;
}

void embed()
{
	unsigned int numSub = msgLen;

	printf("viterbi algorithm :\n");
	int i = 0;
	int j = 0;
	int k = 0;
	for( i = 0 ; i < numSub ; i++ )
	{
		printf("i = %d\n" , i);
		for( j = 0 ; j < w ; j++ )
		{
		printf("	j = %d\n" , j);
			for( k = 0 ; k < numState; k++ )
			{
				printf("		k = %d\n" , k);
				double w0 = ( wght[k] == INFINITY ) ? ( INFINITY ) : ( wght[k] + x[indx] * distortion() );
				double w1 = ( wght[ k ^ matrix[i][j] ] == INFINITY ) ? ( INFINITY ) : ( wght[ k ^ matrix[i][j] ] + (1 - x[indx]) * distortion() );
				path[k][indx] = (unsigned int)( (w1 < w0) ? (1.0) : (0.0) );
				printf("			w0 = %f w1 = %f path[%d][%d] = %d\n",w0 , w1 , k , indx , path[k][indx]);
				newwght[k] = (w1 < w0)?(w1):(w0);
			}
			indx++;
			int l = 0;
			for( l = 0 ; l < numState; l++ )
				wght[l] = newwght[l];
		}

		for( j = 0 ; j < numState; j++ )
		{
			if( j < pow(2 , h - 1) )
				wght[j] = wght[ 2 * j + msg[ indm ] ];
			else
				wght[j] = INFINITY;
		}
		indm++;
	}

	indx--;
	indm--;

	unsigned int state = msg[indm];
	indm--;

	for( i = numSub - 1 ; i >= 0 ; i-- )
	{
		for( j = w - 1 ; j >= 0 ; j-- )
		{
			y[indx] = path[state][indx];
			state = state ^ ( y[indx] * matrix[i][j] );
			indx--;
		}
		if( i == 0 )
			break;
		state = (2 * state + msg[indm])%(int)( pow(2 , h) );

		indm--;
	}

	printf("path matrix :\n");
	for( i = 0 ; i < numState ; i++ )
	{
		for( j = 0 ; j < pathLen ; j++ )
			printf("%d " , path[i][j]);
		printf("\n");
	}

	printf("y :\n");
	for( i = 0 ; i < pathLen ; i++ )
		printf("%d " , y[i] );
	printf("\n");
}

void extract()
{

	printf("extract :\n");
	int i = 0;
	int j = 0;

	unsigned int *out = (unsigned int *)calloc( msgLen , sizeof( unsigned int ) );
	for( i = 0 ; i < msgLen ; i++ )
	{
		unsigned int tmp = 0;
		for( j = 0 ; j < pathLen ; j++ )
			tmp += ( y[j] * Hmatrix[i][j] );	
		out[i] = tmp % 2;
	}


	printf("x :\n");
	for( i = 0 ; i < msgLen ; i++ )
		printf("%d " , out[i]);
	printf("\n");
}

void freePtr()
{
	int i = 0;
	int j = 0;

	free( msg );
	free( x );
	free( y );

	free( wght );
	free( newwght );

	for( i = 0 ; i < h ; i++ )
		free( subMatrix[i] );
	free( subMatrix );

	for( i = 0 ; i < msgLen ; i++ )
	{
		free( matrix[i] );
		free( Hmatrix[i] );
	}
	free( matrix );
	free( Hmatrix );

	for( i = 0 ; i < numState ; i++ )
		free( path[i] );
	free( path );

	return;
}

void main( int argc , char *argv[] )
{
	readInput( argv[1] );

	init();

	embed();

	extract();
	 
	freePtr();
	 
	return;
}
