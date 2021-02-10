classdef POINTSobj
	properties
		coords     = [0,0]
		area       = 0
		Sx         = 0
		Sy         = 0
        panels_in  = []
        panels_out = []
    end
    
    methods 
        function obj = compute_Sx(obj)
            obj.Sx = obj.area * obj.coords(2);
        end
        
        function obj = compute_Sy(obj)
            obj.Sy = obj.area * obj.coords(1);
        end
        
        function obj = change_coords(obj, CG)
            obj.coords = obj.coords - CG;
        end
        
        function cg_p = CG_comp(obj)
           cg_p = [obj.coords(1)*obj.area, obj.coords(2)*obj.area]; 
        end
        
        function obj = rotate(obj,ROT)
            obj.coords = (ROT * obj.coords)';
        end
    end
end