function one_step!(final::Array{ComplexF64},delta_x::Float64,timek::Float64,A::SparseMatrixCSC,B::SparseMatrixCSC,integ)
    rhs = final*[-3,16,-36,48];
    final[1:end,1]=final[1:end,2];
    final[1:end,2]=final[1:end,3];
    final[1:end,3]=final[1:end,4];
    final[1:end,4]=(A+timek*B)\rhs;
    fi =real(sum(final[:,4].*integ));
    fi
end
