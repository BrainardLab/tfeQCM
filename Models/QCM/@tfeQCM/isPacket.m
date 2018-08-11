function packetValidity = isPacket(obj,thePacket)
% packetValidity = isPacket(obj,thePacket)
%
% Function to test if the passed `packet` is well-formed.
%
% 9/14/16   ms      Wrote it.
% 9/21/16   gka     Expanded functionality and reporting

% Check basic validity through parent class
packetValidity = isPacket@tfe(obj,thePacket);

