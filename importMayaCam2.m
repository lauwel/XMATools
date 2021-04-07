function [cam] = importMayaCam2( filename )
% FUNCTION IMPORTMAYACAM2
%   Takes a MayaCam 2.0 file exported by XMALab, parses it, and keeps the
%   pertinant information from the text file:
%     cam.imagesize:
%       a 2-tuple for x and y of the image in pixels
%     cam.M_int: The Intrinsic Camera Matrix
%       focal_x 0 centre_x
%       0 focal_y centre_y
%       0 0 1
%     cam.R: The World to Camera rotation
%     cam.t: the World to Camera translation
%     cam.T_WtC: The World to Camera homogeneous transformation matrix
%       [cam.R, cam.t; 0 0 0 1]
  cam.mayacam2 = filename;

  fid_cam = fopen( cam.mayacam2, 'r' );

  % Import Camera 1
  while 1
    tline = fgetl( fid_cam );
    if ~ischar(tline)
      break;
    elseif regexp(tline, 'image size')
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.imagesize = [str2double(splitline{1,1}), str2double(splitline{2,1})];
    elseif regexp(tline, 'camera matrix')
      cam.M_int = NaN * ones(3,3);
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.M_int(1,:) = [str2double(splitline{1,1}), str2double(splitline{2,1}), str2double(splitline{3,1})];
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.M_int(2,:) = [str2double(splitline{1,1}), str2double(splitline{2,1}), str2double(splitline{3,1})];
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.M_int(3,:) = [str2double(splitline{1,1}), str2double(splitline{2,1}), str2double(splitline{3,1})];
    elseif regexp(tline, 'rotation')
      cam.R = NaN * ones(3,3);
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.R(1,:) = [str2double(splitline{1,1}), str2double(splitline{2,1}), str2double(splitline{3,1})];
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.R(2,:) = [str2double(splitline{1,1}), str2double(splitline{2,1}), str2double(splitline{3,1})];
      tline = fgetl( fid_cam );
      splitline = split(tline, ',');
      cam.R(3,:) = [str2double(splitline{1,1}), str2double(splitline{2,1}), str2double(splitline{3,1})];
    elseif regexp(tline, 'translation')
      cam.t = NaN * ones(3,1);
      tline = fgetl( fid_cam );
      cam.t(1) = str2double(tline);
      tline = fgetl( fid_cam );
      cam.t(2) = str2double(tline);
      tline = fgetl( fid_cam );
      cam.t(3) = str2double(tline);
    end
    clear splitline
  end

  cam.T_WtC = [cam.R, cam.t; 0 0 0 1];

  clear tline

  fclose( fid_cam );

end