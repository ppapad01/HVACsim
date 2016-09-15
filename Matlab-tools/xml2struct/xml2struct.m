function out = xml2struct(xmlfile)
% XML2STRUCT Read an XML file into a MATLAB structure.

% written By Douglas M. Schwarz, douglas.schwarz@kodak.com

xml = xmlread(xmlfile);

children = xml.getChildNodes;
for i = 1:children.getLength
   out(i) = node2struct(children.item(i-1));
end


