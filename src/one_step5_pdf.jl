function one_step5_pdf!(final::Array{ComplexF64},delta_x::Float64,timek::Float64,A::SparseMatrixCSC,B::SparseMatrixCSC,integ)
    rhs = final*[12,-75,200,-300,300];
    final[1:end,1]=final[1:end,2];
    final[1:end,2]=final[1:end,3];
    final[1:end,3]=final[1:end,4];
    final[1:end,4]=final[1:end,5];
    final[1:end,5]=(A+timek*B)\rhs;
    AA=(A-137*I)/(-60*delta_x);
    BB=B/(-60*delta_x);
    ff=(AA+timek*BB)*final[:,5];
    fi=real(sum(ff.*integ));
    fi
end
