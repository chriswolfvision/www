
/***************************************************************************
 The CPU version
 ***************************************************************************/

void cpuFilter(unsigned char *in, unsigned char *out, int rows, int cols)
{

	// General case
	for (int y=1; y<rows-1; ++y)
	for (int x=1; x<cols-1; ++x)
	{
		float f = (
			4.0*in[x*rows+y] +
			2.0*in[(x-1)*rows+y] +
			2.0*in[(x+2)*rows+y] +
			2.0*in[x*rows+y+1] +
			2.0*in[x*rows+y-1] +
			in[(x-1)*rows+y-1] +
			in[(x-1)*rows+y+1] +
			in[(x+1)*rows+y-1] +
			in[(x+1)*rows+y+1]
			)/16.0;
		if (f<0) f=0;
		if (f>255) f=255;
		out[x*rows+y] = (unsigned char) f;
	}
	
	// Borders
	for (int y=0; y<rows; ++y)
	{
		out[0*rows+y] = in[0*rows+y];
		out[(cols-1)*rows+y] = in[(cols-1)*rows+y];
	}
		
	for (int x=0; x<cols; ++x)
	{
		out[x*rows+0] = in[x*rows+0];
		out[x*rows+rows-1] = in[x*rows+rows-1];
	}
}
