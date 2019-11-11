classdef geometry
    properties
        x;
        y;
        r;
        theta;
        resolution;
        circa;
    end
    
    methods
        function obj = geometry(res)
            if (nargin > 0)
                [obj.x,obj.y] = meshgrid(-1:res:1,-1:res:1);
                obj.r=sqrt(obj.x.^2+obj.y.^2);
                obj.theta=atan2(obj.y,obj.x);
                obj.circa = double(gt(1,obj.r));
                obj.circa(obj.circa==0) = NaN;
            end
        end
    end
end

    