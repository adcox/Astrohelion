#include <mgl2/mgl.h>

int main()
    {
        mglData dat(30,40);	// data to for plotting
        for(long i=0;i<30;i++)   for(long j=0;j<40;j++){
            dat.a[i+30*j] = 1/(1+(i-15)*(i-15)/225.+(j-20)*(j-20)/400.);
        }
        mglGraph gr;		// class for plot drawing
        gr.Rotate(50,60);	// rotate axis
        gr.Light(true);		// enable lighting
        gr.Surf(dat);		// plot surface
        gr.Cont(dat,"y");	// plot yellow contour lines
        gr.Axis();			// draw axis
        gr.WriteFrame("sample.png");	// save it
    }