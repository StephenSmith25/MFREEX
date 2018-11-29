#include "connect_segments.h"





int connect_segments(int * segments, int ** connected_segments, int num_segments)
{
	// loop over segments
	int * connect_segments = malloc(num_segments*2*sizeof(int));



	int * segments_copy = malloc(2*num_segments*sizeof(int));
	memcpy(segments_copy,segments,2*num_segments*sizeof(int));

	// set first segment equal to segments[0]
	connect_segments[0]= segments_copy[0];
	connect_segments[1] = segments_copy[1];

	// remove segment from the list
	segments_copy[0] = -1;
	segments_copy[1] = -1;

	for (int i = 1; i < num_segments  ; ++i)
	{
		connect_segments[2*i] = connect_segments[2*(i-1) + 1];
		connect_segments[2*i +1 ] = -1 ;
		int endpoint = connect_segments[2*i];

		// find the endpoint in the list of remaining segments
		for ( int k = 0 ; k < num_segments ; ++k)
		{
			if ( endpoint == segments_copy[2*k] )
			{
				connect_segments[2*i + 1 ] = segments_copy[2*k + 1];
				segments_copy[2*k] = -1;
				segments_copy[2*k+1] = -1;
			}

		}
		if ( connect_segments[2*i+1 ] == -1 )
		{
			for (int k = 0; k < num_segments; ++k)
			{
				if ( endpoint == segments_copy[2*k+1] )
				{
					connect_segments[2*i + 1 ] = segments_copy[2*k];
					segments_copy[2*k] = -1;
					segments_copy[2*k+1] = -1;
				}
			}
		}

	}

	*connected_segments =  connect_segments;

	free(segments_copy);



	return 0;
}
