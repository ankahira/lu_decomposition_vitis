

//------------------------------------------------------------------------------
//
// kernel:  lu
//
// Purpose: Perform LU decomposition of a Matrix
//

#define BUFFER_SIZE 256
#define DATA_SIZE 4096
//TRIPCOUNT identifier
const unsigned int c_len = DATA_SIZE / BUFFER_SIZE;
const unsigned int c_size = BUFFER_SIZE;

/*
    LU Decomposition Kernel Implementation
    Arguments:
    	float * A         --> input Matrix A
    	float * L         --> Output Matrix L
    	float * U         --> Output Matrix U
    	int n			  --> size of matrix A

*/

extern "C"
{
  void krnl_lu(float *A,
               float *L,
               float *U,
               int n

  )
  {
    //perform LU decomposition

    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        if (j < i)
          L[j * n + i] = 0;
        else
        {
          L[j * n + i] = A[j * n + i];
          for (k = 0; k < i; k++)
          {
            L[j * n + i] = L[j * n + i] - L[j * n + k] * U[k * n + i];
          }
        }
      }
      for (j = 0; j < n; j++)
      {
        if (j < i)
          U[i * n + j] = 0;
        else if (j == i)
          U[i * n + j] = 1;
        else
        {
          U[i * n + j] = A[i * n + j] / L[i * n + i];
          for (k = 0; k < i; k++)
          {
            U[i * n + j] = U[i * n + j] - ((L[i * n + k] * U[k * n + j]) / L[i * n + i]);
          }
        }
      }
    }
  }
}


