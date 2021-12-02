# bb

bb is a block sorting (BWT) compressor based on the [bbb compressor by Matt Mahoney](http://mattmahoney.net/dc/#bbb), but without the big block mode.<br>
Each block is BWT transformed and then arithmetic encoded with an indirect context model and 6 SSE stages.<br>
<br>
To Compress: bb.exe c n input output, where n is the block size<br>
To Decompress: bb.exe d input output
