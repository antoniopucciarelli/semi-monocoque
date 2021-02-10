classdef PANELSobj
	properties
		points = zeros(2)
		start  
		finish 
		len    = 0 
		area   = 0 
		flux   = 0 
		flux_p = 0
        link   = []
        thick  = 0
        
    end
    
    methods 
        function obj = get_panels(obj,POINTSarray) 
            obj.start  = POINTSarray(obj.points(1));
            obj.finish = POINTSarray(obj.points(2));
        end
        
        function obj = compute_len(obj)
            obj.len = sqrt((obj.start.coords(1) - obj.finish.coords(1))^2 + ... 
	                       (obj.start.coords(2) - obj.finish.coords(2))^2);
        end
     end
    
end