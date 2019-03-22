/** DynamicArrays.c
 **
 ** Implementation of a (constant growth) Dynamic Array data type                           
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ClassErrors.h"
#include "DynamicArrays.h"


int i = 0;
/***************************************************************************
 CreateDArray -- Allocates memory and initializes DArray state             
        void CreateDArray(DArray *DArrayHead, unsigned int InitialSize)
 Where:
    DArray *DArrayHead       - Pointer to a storage block which holds the
                               data structure entry 
    unsigned int InitialSize - The initial data block size, could be 0
    returns: void            - nothing is returned.
    errors:                  - This routine will print an error message to 
                               stderror and exit with an error code
 **************************************************************************/
void CreateDArray(DArray *DArrayHead, unsigned int InitialSize) {
   
  /*Initializing array*/
  DArrayHead->EntriesUsed = 0;

  /*Initializing capacity to u_int InitialSize */
  DArrayHead->Capacity = InitialSize;
  
  /*Allocate storage space for payload according to struc Data*/
  if (InitialSize==0)
  {
  	DArrayHead->Payload = malloc(0);
  }
  
  else 
  {
	DArrayHead->Payload = malloc(InitialSize * sizeof(Data));
  }

  /*If memory not allocated print error to stderr and exit with a code */
  if (DArrayHead->Payload == NULL)
  {
  	printf("ERROR: malloc failed, DynamicArrays.c Line:55 \n");
        exit(MALLOC_ERROR);
  }
 
} /* CreateDArray() */


/***************************************************************************
  PushToDArray -- Adds data at the end of the dynamic array and updates its 
                  state in the array header.  If full, reallocates memory 
                  space according to growth policy                          
        unsigned int PushToDArray(DArray *DArrayHead, Data *PayloadPtr)     
Where:
    DArray *DArrayHead      - Pointer to a storage block which holds the
                              data structure entry 
    Data * PayloadPtr       - 
    returns: unsigned int   - The index of the last element inserted
    errors:                 - This routine will print an error message to 
                              stderror and exit with an error code 
 **************************************************************************/
unsigned int PushToDArray(DArray *DArrayHead, Data *PayloadPtr)
{

  /*Check array size */
  if (DArrayHead->EntriesUsed >= DArrayHead->Capacity)
  {
   
   /*If required, increase array  size */
   DArrayHead->Capacity += GROWTH_AMOUNT;
   
   /*Re-allocate storage for elements */
   DArrayHead->Payload = realloc(DArrayHead->Payload, sizeof(Data) * (DArrayHead->Capacity));
   
   /*Print error message if memeory is not reallocated */
   if (DArrayHead->Payload == NULL) 
   {
   	 printf("ERROR: realloc failed, DynamicArrays.c Line:91 \n");
         exit(REALLOC_ERROR);
   }
  }
  
  /*Copy data from PayloadPtr to DArrayHead */
  memcpy(&(DArrayHead->Payload[DArrayHead->EntriesUsed]), PayloadPtr, sizeof(Data));
  
  /*Increase number of elements*/
  (DArrayHead->EntriesUsed)++;

  /*Return array index of last element inserted */
  return (DArrayHead->EntriesUsed - 1);

} /* PushToDArray() */


/*************************************************************************
 DestroyDArray -- Returns the memory back to the heap and updates the   
                  state in the array header      
  void DestroyDArray(DArray *DArrayHead)
  Where:
    DArray *DArrayHead  - Pointer to a storage block which holds the
                          data structure entry 
    returns: void       - nothing is returned.
    errors:             - No errors are reported
*************************************************************************/
void DestroyDArray(DArray *DArrayHead)
{
  /*DArrayHead -> Number of entries used = 0 */
  DArrayHead->EntriesUsed = 0;
  
  /*DArrayHead -> Capacity = 0 */
  DArrayHead->Capacity = 0;

  /*De-allocate storage */
  free(DArrayHead->Payload);
  DArrayHead->Payload = NULL;
    
} /* DestroyDArray() */
