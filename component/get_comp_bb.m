function bb = get_comp_bb(img, max_comp)
% GET_COMP_BB  Determine the bounding box of each component in the labelled img
%
% BB = get_comp_bb(IMG, MAX_COMP)
%
% This function takes a labelled component image as input, and determines the 
% left, top, right, and bottom pixel co-ordinates (bounding box), of each 
% component.
%
% IMG should be 2D component label matrix, with each pixel of a component is
% assigned the same number, and the numbers are assumed to be contiguous from 
% 1...MAX_COMP (see bwlabel)
%
% MAX_COMP should be the maximum component number in the image.  It is passed
% instead of calculated directly for efficiency reasons (should be available at
% call time if bwlabel has been used to generate the matrix).
%
% BB will be an Nx4 matrix where each row specifies the left,top,right, and
% bottom pixel location for a single component
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: get_comp_bb.m,v 1.1 2006-09-20 21:45:33 scottl Exp $
%
% REVISION HISTORY
% $Log: get_comp_bb.m,v $
% Revision 1.1  2006-09-20 21:45:33  scottl
% initial check-in.
%

% LOCAL VARS %
%%%%%%%%%%%%%%
classname = 'uint16';  %what type of values should bb take on.


% CODE START %
%%%%%%%%%%%%%%
if nargin ~= 2
    error('incorrect number of arguments specified!');
end

bb = zeros(max_comp, 4, classname);

%fill in the L T R B pixel boundary values for each component
[R, C] = find(img ~= 0);
for ii=1:length(R)
    loc = img(R(ii), C(ii));
    %is this the first time seeing this component?
    if bb(loc,1) == 0
        %set all 4 positions for this component to this location
        bb(loc,1) = C(ii); bb(loc,2) = R(ii); 
        bb(loc,3) = C(ii); bb(loc,4) = R(ii);
    else
        %update the R position and see if we should update T or B
        bb(loc,3) = C(ii);
        if R(ii) < bb(loc,2)
            bb(loc,2) = R(ii);
        end
        if R(ii) > bb(loc,4)
            bb(loc,4) = R(ii);
        end
    end
end
