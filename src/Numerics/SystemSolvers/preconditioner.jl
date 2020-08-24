export AbstractPreconditioner, ColumnwiseLUPreconditioner, preconditioner_update!, preconditioner_solve!


abstract type AbstractPreconditioner end



mutable struct ColumnwiseLUPreconditioner{AT} <: AbstractPreconditioner 
    A::DGColumnBandedMatrix
    Q::AT
    PQ::AT
    counter::Int64
    update_freq::Int64
end

# function ColumnwiseLUPreconditioner(op, dg, Q0, args...)
#     single_column = false
#     Q = similar(Q0)
#     PQ= similar(Q0)

#     # TODO: can we get away with just passing the grid?
#     A = banded_matrix(
#         op,
#         dg,
#         Q,
#         PQ,
#         args...;
#         single_column = single_column,
#     )

#     band_lu!(A)
#     ColumnwiseLUPreconditioner(A, Q, PQ, 0)
# end


function ColumnwiseLUPreconditioner(dg, Q0, update_freq=100)
    single_column = false
    Q = similar(Q0)
    PQ = similar(Q0)

    # TODO: can we get away with just passing the grid?
    A = banded_matrix(
        dg,
        Q;
        single_column = single_column,
    )

    band_lu!(A)
    ColumnwiseLUPreconditioner(A, Q, PQ, -1, update_freq)
end

"""
update the DGColumnBandedMatrix
"""
function preconditioner_update!(op, dg, preconditioner::ColumnwiseLUPreconditioner, args...)
    
    if preconditioner.counter >= 0 && (preconditioner.counter < preconditioner.update_freq)
        return
    end

    A = preconditioner.A
    Q = preconditioner.Q
    PQ = preconditioner.PQ

    # TODO: can we get away with just passing the grid?
    update_banded_matrix!(
        A,
        op,
        dg,
        Q,
        PQ,
        args...
    )

    band_lu!(A)
    preconditioner.counter = 0
end

"""
Inplace solve Q = Pinv * Q
"""
function preconditioner_solve!(preconditioner::ColumnwiseLUPreconditioner, Q)
    A = preconditioner.A
    band_forward!(Q, A)
    band_back!(Q, A)


    preconditioner.counter += 1

end


