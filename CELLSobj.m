classdef CELLSobj
	properties
		coeffs
        panels = [];
		PANELS  
	end
	
	methods 
		function get_panels(obj,panel_array)
			for ii = 1:length(obj.panels)
				obj.PANELS(ii) = panel_array(obj.panels(ii));
			end
		end
	end

end