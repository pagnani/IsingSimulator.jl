using PrecompileTools

@setup_workload begin
    @compile_workload begin
        redirect_stdout(devnull) do
            N=10
            ns=10;
            isingview((N, N), βc2d - 0.01, nsweep=ns, fps=10^8, x0=ones(Int, N^2))
        end
    end
end