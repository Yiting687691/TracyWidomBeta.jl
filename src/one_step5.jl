function one_step5!(final::Array{ComplexF64},delta_x::Float64,timek::Float64,A::SparseMatrixCSC,B::SparseMatrixCSC,integ)
    rhs = final*[12,-75,200,-300,300];
    final[1:end,1]=final[1:end,2];
    final[1:end,2]=final[1:end,3];
    final[1:end,3]=final[1:end,4];
    final[1:end,4]=final[1:end,5];
    final[1:end,5]=(A+timek*B)\rhs;
    fi =real(sum(final[:,5].*integ));
    fi
end
