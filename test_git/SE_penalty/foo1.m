function a = foo1(varargin)
mear = varargin(1:end);
    a = foo2(mear{:});
end

function a = foo2(varargin)
varargin
a = feval(foo,varargin{:});
end

function a = foo(b,c)
a = b + c;
end