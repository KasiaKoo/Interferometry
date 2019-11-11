classdef PSI
    properties
        al;
        t;
        r;
        wrap;
        phase;
        c_nm;
    end
    
    methods
        function obj = PSI(test,reference,algorithm,g,how)
            if (nargin>0)
                obj.t = test;
                obj.r = reference;
                obj.al = algorithm;
                if algorithm == 'example'
                    obj=obj.example(g,how);
                end
                
            end               
           
        end
        
        function obj = example(obj,g,how)
            I1=abs(obj.r+obj.t).^2;
            I2=abs(exp(1i*2*pi/3).*obj.r + obj.t).^2;
            I3=abs(exp(1i*4*pi/3).*obj.r + obj.t).^2;
            obj.wrap=g.circa.*atan2(sqrt(3)/2*(I2-I3),(I1-(I2+I3)/2));
            obj = obj.unwrap_phase(how);
        end
            
        function obj = unwrap_phase(obj, how,g)
            if how == 'centre'
                q=obj.wrap;
                q(200:end,:)=unwrap(q(200:end,:),[],1);
                q=flipud(q);
                q(200:end,:)=unwrap(q(200:end,:),[],1);
                q(:,200:end)=unwrap(q(:,200:end),[],2);
                q=fliplr(q);
                q(:,200:end)=unwrap(q(:,200:end),[],2);
                obj.phase = rot90(q,2);
                obj = obj.ref_takeaway(g);
                obj = obj.zernike_coefficients(g);
            elseif how =='edge'
                obj.phase = unwrap(unwrap(p,[],2),[],1);
                obj = obj.ref_takeaway(g);
                obj = obj.zernike_coefficients(g);
            end
        end
        
        function obj = ref_takeaway(obj,g)
            obj.phase = g.circa.*(obj.phase + obj.r.height_error);
                 
        end
        
        function obj = zernike_coefficients(obj,g)
            Z_00 = 1.*g.circa;
            Z_11 = 2.*g.r.*cos(g.theta).*g.circa;
            Z_1n1 = 2.*g.r.*sin(g.theta).*g.circa;
            Z_20 = sqrt(3).*(2.*(g.r.^2) - 1).*g.circa;
            Z_22 = sqrt(6).*(g.r.^2).*cos(2.*g.theta).*g.circa;
            Z_2n2 = sqrt(6).*(g.r.^2).*sin(2.*g.theta).*g.circa;
            Z_31 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*cos(g.theta).*g.circa;
            Z_3n1 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*sin(g.theta).*g.circa;
            Z_33 = 2*sqrt(2).*(g.r.^3).*cos(3.*g.theta).*g.circa;
            Z_3n3 = 2*sqrt(2).*(g.r.^3).*sin(3.*g.theta).*g.circa;
            Z_40 = sqrt(5).*(6.*(g.r.^4) - 6.*(g.r.^2)+1).*g.circa;
            c_lst = linspace(0,0,11);
            abb_options = {Z_00, Z_11, Z_1n1, Z_20, Z_22, Z_2n2, Z_31, Z_3n1, Z_33, Z_3n3, Z_40};
            for i=1:11
                abe = abb_options(i);
                c_lst(i) = sum(obj.phase(:).*abe(:))/sum(g.circa(:));
            end
            obj.c_nm = c_lst;
        end    
    end
    
end

    




