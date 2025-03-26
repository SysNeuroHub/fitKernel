%for figure 1

addpath(genpath('/home/daisuke/Documents/git/neurostim'));
addpath(genpath('/home/daisuke/Documents/git/marmolab-stimuli'));

% /home/daisuke/Documents/git/marmolab-common/facecal.m
  rigprefs = getpref('marmolab','rigprefs');

tmp = load(rigprefs.facelib);
tmp = struct2cell(tmp);
tmp = tmp([7,10,13,17:20,24,25,27]); % these faces seem most centered

c = marmolab.rigcfg('debug',true);
faces = neurostim.stimuli.texture(c,'faces');

N = length(tmp);

for id = 1:N
  img = tmp{id};

  % gaussian window...
  sz = size(img);
  img(:,:,4) = uint8(255*faces.mkwin(max(sz(:)),0.15)); % alpha channel: 0 = transparent, 255 = opaque

  faces.add(id,img);

  subplot(4,3,id);
  image(img(:,:,1:3),'AlphaData',img(:,:,4));
  axis off equal
end


thisid = 2;
f = figure('position', [0 0 200 200]);
img = tmp{thisid};
sz = size(img);
img(:,:,4) = uint8(255*faces.mkwin(max(sz(:)),0.15)); % alpha channel: 0 = transparent, 255 = opaque
image(img(:,:,1:3),'AlphaData',img(:,:,4));
axis off equal
screen2png('MarmoFace');

