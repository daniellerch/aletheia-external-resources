default:
	mkoctfile --mex jpeg_read.c -ljpeg -fno-lto
	mkoctfile --mex jpeg_write.c -ljpeg -fno-lto
	cp -f *.mex ../octave/
clean:
	rm -f *.o *.mex
