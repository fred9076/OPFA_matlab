function varargout = geotiffread_modified(filename, varargin)
%GEOTIFFREAD Read GeoTIFF file
%
%   [A, R] = GEOTIFFREAD(FILENAME) reads a georeferenced grayscale, RGB, or
%   multispectral image or data grid from the GeoTIFF file specified by 
%   FILENAME into A and constructs a spatial referencing object, R.
%
%   [X, CMAP, R] = GEOTIFFREAD(FILENAME) reads an indexed image into X and
%   the associated colormap into CMAP, and constructs a spatial referencing
%   object, R. Colormap values in the image file are rescaled into the
%   range [0,1].
%
%   FILENAME is a string scalar or character vector that specifies the name
%   of the GeoTIFF file. FILENAME can include the folder name. Otherwise,
%   the file must be in the current folder or in a folder on the MATLAB
%   path. If the named file includes the extension '.TIF' or '.TIFF'
%   (either upper or lower case), you can omit the extension from FILENAME.
%
%   A is a two-dimensional array if the file contains a grayscale image or
%   data grid. A is an M-by-N-by-P array if the file contains a color
%   image, multispectral image, hyperspectral image, or data grid. The
%   class of A depends on the storage class of the pixel data in the file
%   which is related to the BitsPerSample property as returned by the
%   IMFINFO function.
%
%   R is a geographic raster reference object if the image or data grid is
%   referenced to a geographic coordinate system, or a map raster reference
%   object if it is referenced to a projected coordinate system.
%
%   [A, REFMAT, BBOX] = GEOTIFFREAD(FILENAME) reads a georeferenced
%   grayscale, RGB, or multispectral image or data grid into A; the
%   corresponding referencing matrix into REFMAT; and the bounding box into
%   BBOX.
%
%   [X, CMAP, REFMAT, BBOX] = GEOTIFFREAD(FILENAME) reads an indexed image
%   into X, the associated colormap into CMAP, the referencing matrix into
%   REFMAT, and the bounding box into BBOX.  The referencing  matrix must
%   be unambiguously defined by the GeoTIFF file, otherwise  it and the
%   bounding box are returned empty.
%
%   [...] = GEOTIFFREAD(FILENAME, IDX) reads one image from a multi-image
%   GeoTIFF file. IDX is an integer value that specifies the order that the
%   image appears in the file. For example, if IDX is 3, GEOTIFFREAD reads
%   the third image in the file. If you omit this argument, GEOTIFFREAD
%   reads the first image in the file.
%
%   [...] = GEOTIFFREAD(URL, ...) reads the GeoTIFF image from a URL. The
%   URL must include the protocol type (e.g., "http://").
%
%   Note
%   ----
%   GEOTIFFREAD imports pixel data using the TIFF-reading capabilities of
%   the MATLAB function IMREAD and likewise shares any limitations of
%   IMREAD.  Consult the IMREAD documentation for specific information on
%   TIFF image support.
%
%   Example
%   -------
%   % Read and display the Boston GeoTIFF image. 
%   % Includes material (c) GeoEye, all rights reserved.
%   [boston, R] = geotiffread('boston.tif');
%   figure
%   mapshow(boston, R)
%   axis image off
%   
%   See also GEOSHOW, GEOTIFFINFO, GEOTIFFWRITE, IMREAD, MAPSHOW

% Copyright 1996-2017 The MathWorks, Inc.

% Verify the input and output argument count.
narginchk(1,2);
nargoutchk(0,4);

% Parse the inputs.
filename = convertStringsToChars(filename);
[filename, url, idx] = parseInputs(filename, varargin);

% Read the info fields from the filename.
info = geotiffinfo(filename);

% Read the image from the filename.
[A, cmap] = imread(filename,idx);

% If the RefMatrix is empty, try to obtain the spatial information from a
% corresponding worldfile.
if isempty(info.RefMatrix)
    worldfilename = getworldfilename(filename);
    info = getSpatialInfoFromWorldfile(worldfilename, info, size(A));
end

% Delete temporary file from Internet download.
if (url)
    deleteDownload(filename);
end

% Assign output arguments.
varargout = assignOutputArguments(A, cmap, info, nargout);

%--------------------------------------------------------------------------

function [filename, url, idx] = parseInputs(filename, inputs)
% Parse the inputs from the cell array, INPUTS.

% Verify the filename and obtain the full pathname.
extensions = {'tif', 'tiff'};
[filename, url] = internal.map.checkfilename(filename, extensions, mfilename, 1, true);

% Check and set the image index number.
if ~isempty(inputs)
  idx = inputs{1};
  attributes = {'real' 'scalar' 'positive'};
  validateattributes(idx, {'numeric'}, attributes, mfilename, 'IDX', 2);
else
  idx = 1;
end

%--------------------------------------------------------------------------

function info = getSpatialInfoFromWorldfile(worldfilename, info, rasterSize)
% Obtain the referencing matrix from a world file, if it exits. If so,
% compute the bounding box and construct a spatial referencing object from
% the referencing matrix and update the fields of the INFO structure.
% WORLDFILENAME is a string denoting the name of the world file.

if exist(worldfilename,'file')
    % Obtain the referencing matrix from the world file.
    refmat = worldfileread(worldfilename);
    
    % Calculate the spatial referencing object and bounding box from the
    % referencing matrix if it is not empty.
    if ~isempty(refmat)               
        if strcmp(info.ModelType, 'ModelTypeGeographic')
            R = refmatToGeoRasterReference(refmat, rasterSize);
        elseif strcmp(info.ModelType, 'ModelTypeProjected')
            R = refmatToMapRasterReference(refmat, rasterSize);
        else
            R = [];
        end
        info.BoundingBox = mapbbox(refmat, rasterSize);
        info.RefMatrix = refmat;
        info.SpatialRef = R;
    end   
end

%--------------------------------------------------------------------------

function outputs = assignOutputArguments(A, cmap, info, numOutputs)
% Assign the output arguments based on the number of arguments requested.

outputs{1} = A;
switch numOutputs
    case 2
        if strcmp(info.ColorType, 'indexed')
            % [X, CMAP] = GEOTIFFREAD(...)
            outputs{2} = cmap;
        else
            % [A, R] = GEOTIFFREAD(...)
            outputs{2} = info.SpatialRef;
        end
        
    case 3
        if strcmp(info.ColorType, 'indexed')
            % [X, CMAP, R] = GEOTIFFREAD(...)
            outputs{2} = cmap;
            outputs{3} = info.SpatialRef;
        else
            % [A, REFMAT, BBOX] = GEOTIFFREAD(...)          
            outputs{2} = info.RefMatrix;
            outputs{3} = info.BoundingBox;
        end
        
    case 4
        % [X, CMAP, REFMAT, BBOX] = GEOTIFFREAD(...)
        outputs{2} = cmap;
        outputs{3} = info.RefMatrix;
        outputs{4} = info.BoundingBox;
end
