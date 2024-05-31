function data = readRays(fname,varargin)

if nargin < 2
    nstart = 1;
    nframes = Inf;
else
    nstart=varargin{1};
    nframes=varargin{2};
end


fid = fopen(fname);
stop = 0;
nread = 0;
kframe = 0;

while and((nread < nframes),(~feof(fid)))
    kframe=kframe+1;
    nrays = fread(fid,1,'uint32');
    if kframe >= nstart
        kframe
        nrays
        nread=nread+1
        for kray=1:nrays
            camID(kray) = fread(fid,1,'uint8');
            rayID(kray) = fread(fid,1,'uint16');
            P(kray,1:3) = fread(fid,3,'float32');
            V(kray,1:3) = fread(fid,3,'float32');
        end
        if nrays~=0
            data(nread).camID=camID;
            data(nread).rayID=rayID;
            data(nread).P=P;
            data(nread).V=V;
        else
            data(nread).camID=[];
            data(nread).rayID=[];
            data(nread).P=[];
            data(nread).V=[];
        end

        clear camID rayID P V;

    else
        fseek(fid,27*nrays,0);
    end
end

fclose(fid);