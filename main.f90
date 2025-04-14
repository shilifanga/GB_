 program mass_matrix_assembly
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: i, j, k, ios
  real(dp) :: a, b, h, threshold, err, temp, sum_val, xq
  integer, parameter :: numQuad = 5
  real(dp), dimension(numQuad) :: quad_points, quad_weights
  integer, parameter :: N = 10
  integer, parameter :: nM = N + 5
  real(dp) :: M(Nm,Nm)
  real(dp), dimension(15,15) :: xia  ! 假定 xiadate.txt 文件为 13x15 数字数据
  integer :: row, col
  integer :: cols
 
  ! 设置问题参数
  a = 0.0_dp
  b = 1.0_dp

  h = (b - a) / real(N, dp)

  ! 积分点与权重（高斯积分节点）
  quad_points = (/ -0.9061798459386640_dp, -0.5384693101056830_dp, 0.0_dp, 0.5384693101056830_dp, 0.9061798459386640_dp /)
  quad_weights = (/  0.2369268850561890_dp,  0.4786286704993660_dp, 0.5688888888888889_dp, 0.4786286704993660_dp, 0.2369268850561890_dp /)

  ! 分配矩阵 M (尺寸为 (N+5) x (N+5) )
   
  M = 0.0_dp

  ! 组装矩阵，包含内部点及边界处理
  call assmblep4(N, h, quad_points, quad_weights, M)

  ! Meanwhile, the ghost cells are also included in the assembly of the matrix.
  M(1,1) = 1_dp
  M(2,2) = 1_dp
  M(nM-1,nM-1) = 1_dp
  M(nM,nM) = 1_dp
  ! 读入参考数据文件 xiadate.txt （假定为 15 行 15 列数据）
  !open(newunit=i, file='xiadate.txt', status='old', action='read', iostat=ios)
  !if (ios /= 0) then
  !  print *, 'Error opening file xiadate.txt'
  !  stop
  !end if
  !do row = 1, 15
  !   read(i,*, iostat=ios) (xia(row, col), col=1,15)
  !   if (ios /= 0) then
  !     print *, 'Error reading file xiadate.txt at row ', row
  !     stop
  !   end if
  !end do
  !close(i)

  ! 检查计算结果与文件数据的误差（仅对第 3 到 13 行进行比较）
  !threshold = 1.0e-15_dp
  !do k = 1, 15
  !   err = 0.0_dp
  !   do j = 1, 15
  !      err = max(err, abs( h * xia(k,j) - M(k,j) ))
  !   end do
  !   if (err > threshold) then
  !      print *, 'mass matrix is not correct! Row ', k, ' error = ', err
  !   end if
  !end do

    print *, "Matrix:"
    do i = 1, nM
        do j = 1, nM
            write(*, '(F10.8)', advance='no') M(i, j)
        end do
        print *
    end do
  print *, 'Mass matrix assembly complete.'
  
contains

  !===========================================================================
  ! 子程序：assmblep4
  ! 功能：组装质量矩阵 M，尺寸为 (N+5) x (N+5)
  !===========================================================================
  subroutine assmblep4(N, h, quad_points, quad_weights, M)
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in) :: h
    real(dp), intent(in), dimension(:) :: quad_points, quad_weights
    real(dp), intent(inout), dimension(N+5, N+5) :: M

    integer, parameter :: numQuad = 5
    integer :: n1, q
    real(dp) :: sum_val, xq

    ! 注意：以下这些变量计算积分结果，实际上与 n1 无关（在均匀网格下为常数）
    real(dp) :: N1_N1M3_1, N1_N1M3_2
    real(dp) :: N1_N1M2_1, N1_N1M2_2, N1_N1M2_3
    real(dp) :: N1_N1M1_1, N1_N1M1_2, N1_N1M1_3, N1_N1M1_4
    real(dp) :: N1_N1_1, N1_N1_2, N1_N1_3, N1_N1_4, N1_N1_5
    real(dp) :: N1_N1P1_1, N1_N1P1_2, N1_N1P1_3, N1_N1P1_4
    real(dp) :: N1_N1P2_1, N1_N1P2_2, N1_N1P2_3
    real(dp) :: N1_N1P3_1, N1_N1P3_2

    !---------------------------
    ! STEP1：内部点（n1 = 6 ... N）
    !---------------------------
    do n1 = 6, N

      ! M(n1, n1-4)
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 4.0_dp
        ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift8(5)
        sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(5, xq + 8.0_dp)
      end do
      M(n1, n1-4) = (h/2.0_dp) * sum_val

      ! M(n1, n1-3) = N1_N1M3_1 + N1_N1M3_2
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 4.0_dp
        ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift6(4)
        sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(4, xq + 6.0_dp)
      end do
      N1_N1M3_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 2.0_dp
        ! 使用 GD_basis_p4(2) 与 GD_basis_p4leftshift6(5)
        sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(5, xq + 6.0_dp)
      end do
      N1_N1M3_2 = (h/2.0_dp) * sum_val
      M(n1, n1-3) = N1_N1M3_1 + N1_N1M3_2

      ! M(n1, n1-2) = N1_N1M2_1 + N1_N1M2_2 + N1_N1M2_3
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 4.0_dp
        ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift4(3)
        sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(3, xq + 4.0_dp)
      end do
      N1_N1M2_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 2.0_dp
        ! 使用 GD_basis_p4(2) 与 GD_basis_p4leftshift4(4)
        sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(4, xq + 4.0_dp)
      end do
      N1_N1M2_2 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q)
        ! 使用 GD_basis_p4(3) 与 GD_basis_p4leftshift4(5)
        sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(5, xq + 4.0_dp)
      end do
      N1_N1M2_3 = (h/2.0_dp) * sum_val
      M(n1, n1-2) = N1_N1M2_1 + N1_N1M2_2 + N1_N1M2_3

      ! M(n1, n1-1) = N1_N1M1_1 + N1_N1M1_2 + N1_N1M1_3 + N1_N1M1_4
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 4.0_dp
        ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift2(2)
        sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(2, xq + 2.0_dp)
      end do
      N1_N1M1_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 2.0_dp
        ! 使用 GD_basis_p4(2) 与 GD_basis_p4leftshift2(3)
        sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(3, xq + 2.0_dp)
      end do
      N1_N1M1_2 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q)
        ! 使用 GD_basis_p4(3) 与 GD_basis_p4leftshift2(4)
        sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(4, xq + 2.0_dp)
      end do
      N1_N1M1_3 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 2.0_dp
        ! 使用 GD_basis_p4(4) 与 GD_basis_p4leftshift2(5)
        sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(5, xq + 2.0_dp)
      end do
      N1_N1M1_4 = (h/2.0_dp) * sum_val
      M(n1, n1-1) = N1_N1M1_1 + N1_N1M1_2 + N1_N1M1_3 + N1_N1M1_4

      ! M(n1, n1) = N1_N1_1 + N1_N1_2 + N1_N1_3 + N1_N1_4 + N1_N1_5
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 4.0_dp
        ! 使用 GD_basis_p4(1) 与 GD_basis_p4(1)
        sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(1, xq)
      end do
      N1_N1_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 2.0_dp
        sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(2, xq)
      end do
      N1_N1_2 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q)
        sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(3, xq)
      end do
      N1_N1_3 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 2.0_dp
        sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(4, xq)
      end do
      N1_N1_4 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 4.0_dp
        sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(5, xq)
      end do
      N1_N1_5 = (h/2.0_dp) * sum_val
      M(n1, n1) = N1_N1_1 + N1_N1_2 + N1_N1_3 + N1_N1_4 + N1_N1_5

      ! M(n1, n1+1) = N1_N1P1_1 + N1_N1P1_2 + N1_N1P1_3 + N1_N1P1_4
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) - 2.0_dp
        ! 使用 GD_basis_p4(2) 与 GD_basis_p4rightshift2(1)
        sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(1, xq - 2.0_dp)
      end do
      N1_N1P1_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q)
        ! 使用 GD_basis_p4(3) 与 GD_basis_p4rightshift2(2)
        sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(2, xq -2.0_dp)
      end do
      N1_N1P1_2 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 2.0_dp
        ! 使用 GD_basis_p4(4) 与 GD_basis_p4rightshift2(3)
        sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(3, xq - 2.0_dp)
      end do
      N1_N1P1_3 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 4.0_dp
        ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift2(4)
        sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(4, xq - 2.0_dp)
      end do
      N1_N1P1_4 = (h/2.0_dp) * sum_val
      M(n1, n1+1) = N1_N1P1_1 + N1_N1P1_2 + N1_N1P1_3 + N1_N1P1_4

      ! M(n1, n1+2) = N1_N1P2_1 + N1_N1P2_2 + N1_N1P2_3
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q)
        ! 使用 GD_basis_p4(3) 与 GD_basis_p4rightshift4(1)
        sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(1, xq - 4.0_dp)
      end do
      N1_N1P2_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 2.0_dp
        ! 使用 GD_basis_p4(4) 与 GD_basis_p4rightshift4(2)
        sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(2, xq - 4.0_dp)
      end do
      N1_N1P2_2 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 4.0_dp
        ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift4(3)
        sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(3, xq - 4.0_dp)
      end do
      N1_N1P2_3 = (h/2.0_dp) * sum_val
      M(n1, n1+2) = N1_N1P2_1 + N1_N1P2_2 + N1_N1P2_3

      ! M(n1, n1+3) = N1_N1P3_1 + N1_N1P3_2
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 2.0_dp
        ! 使用 GD_basis_p4(4) 与 GD_basis_p4rightshift6(1)
        sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(1, xq - 6.0_dp)
      end do
      N1_N1P3_1 = (h/2.0_dp) * sum_val
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 4.0_dp
        ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift6(2)
        sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(2, xq - 6.0_dp)
      end do
      N1_N1P3_2 = (h/2.0_dp) * sum_val
      M(n1, n1+3) = N1_N1P3_1 + N1_N1P3_2

      ! M(n1, n1+4)
      sum_val = 0.0_dp
      do q = 1, numQuad
        xq = quad_points(q) + 4.0_dp
        ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift8(1)
        sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(1, xq - 8.0_dp)
      end do
      M(n1, n1+4) = (h/2.0_dp) * sum_val

    end do  ! 内部点循环

    !---------------------------
    ! STEP2：边界点处理
    !---------------------------
    ! 处理第5行
    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 7.0_dp)
      ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift8(5)
      sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(5, xq + 8.0_dp)
    end do
    M(5,1) = (h/2.0_dp * sum_val) / 2.0_dp

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 7.0_dp)
      ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift6(4)
      sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(4, xq + 6.0_dp)
    end do
    M(5,2) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1M3_2

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 7.0_dp)
      ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift4(3)
      sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(3, xq + 4.0_dp)
    end do
    M(5,3) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1M2_2 + N1_N1M2_3

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 7.0_dp)
      ! 使用 GD_basis_p4(1) 与 GD_basis_p4leftshift2(2)
      sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(2, xq + 2.0_dp)
    end do
    M(5,4) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1M1_2 + N1_N1M1_3 + N1_N1M1_4

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 7.0_dp)
      ! 使用 GD_basis_p4(1) 与 GD_basis_p4(1)
      sum_val = sum_val + quad_weights(q) * phi_p4(1, xq) * phi_p4(1, xq)
    end do
    M(5,5) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1_2 + N1_N1_3 + N1_N1_4 + N1_N1_5

    ! 对于 test function 的其他赋值，直接复制内部矩阵部分
    M(5,6) = M(6,7)
    M(5,7) = M(6,8)
    M(5,8) = M(6,9)
    M(5,9) = M(6,10)

    ! 处理第4行（对应 φ1）
    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 3.0_dp)
      ! 使用 GD_basis_p4(2) 与 GD_basis_p4leftshift6(5)
      sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(5, xq + 6.0_dp)
    end do
    M(4,1) = (h/2.0_dp * sum_val) / 2.0_dp

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 3.0_dp)
      ! 使用 GD_basis_p4(2) 与 GD_basis_p4leftshift4(4)
      sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(4, xq + 4.0_dp)
    end do
    M(4,2) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1M2_3

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 3.0_dp)
      ! 使用 GD_basis_p4(2) 与 GD_basis_p4leftshift2(3)
      sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(3, xq + 2.0_dp)
    end do
    M(4,3) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1M1_3 + N1_N1M1_4

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 3.0_dp)
      ! 使用 GD_basis_p4(2) 与 GD_basis_p4(2)
      sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(2, xq)
    end do
    ! 修改处：将 N1_N1M1_3,N1_N1M1_4,N1_N1M1_5 替换为 N1_N1_3,N1_N1_4,N1_N1_5
    M(4,4) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1_3 + N1_N1_4 + N1_N1_5

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 3.0_dp)
      ! 使用 GD_basis_p4(2) 与 GD_basis_p4rightshift2(1)
      sum_val = sum_val + quad_weights(q) * phi_p4(2, xq) * phi_p4(1, xq - 2.0_dp)
    end do
    M(4,5) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1P1_2 + N1_N1P1_3 + N1_N1P1_4

    M(4,6) = M(5,7)
    M(4,7) = M(5,8)
    M(4,8) = M(5,9)

    ! 处理第3行（对应 φ0）
    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4leftshift4(5)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(5, xq + 4.0_dp)
    end do
    M(3,1) = (h/2.0_dp * sum_val) / 2.0_dp

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4leftshift2(4)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(4, xq + 2.0_dp)
    end do
    M(3,2) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1M1_4

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4(3)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(3, xq)
    end do
    M(3,3) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1_4 + N1_N1_5

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4rightshift2(2)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(2, xq - 2.0_dp)
    end do
    M(3,4) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1P1_3 + N1_N1P1_4

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4rightshift4(1)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(1, xq - 4.0_dp)
    end do
    M(3,5) = (h/2.0_dp * sum_val) / 2.0_dp + N1_N1P2_2 + N1_N1P2_3

    M(3,6) = M(4,7)
    M(3,7) = M(4,8)

    ! 右边界的处理
    M(N+1, N-3) = M(N, N-4)
    M(N+1, N-2) = M(N, N-3)
    M(N+1, N-1) = M(N, N-2)
    M(N+1, N)   = M(N, N-1)

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 7.0_dp)
      ! 使用 GD_basis_p4(5) 与 GD_basis_p4(5)
      sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(5, xq)
    end do
    M(N+1, N+1) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1_1 + N1_N1_2 + N1_N1_3 + N1_N1_4

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 7.0_dp)
      ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift2(4)
      sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(4, xq - 2.0_dp)
    end do
    M(N+1, N+2) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1P1_1 + N1_N1P1_2 + N1_N1P1_3

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 7.0_dp)
      ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift4(3)
      sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(3, xq - 4.0_dp)
    end do
    M(N+1, N+3) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1P2_1 + N1_N1P2_2

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 7.0_dp)
      ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift6(2)
      sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(2, xq - 6.0_dp)
    end do
    M(N+1, N+4) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1P3_1

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 7.0_dp)
      ! 使用 GD_basis_p4(5) 与 GD_basis_p4rightshift8(1)
      sum_val = sum_val + quad_weights(q) * phi_p4(5, xq) * phi_p4(1, xq - 8.0_dp)
    end do
    M(N+1, N+5) = (h/2.0_dp * sum_val)/2.0_dp

    ! 对于 M(N+2, *)（对应 φ9）
    M(N+2, N-2) = M(N, N-4)
    M(N+2, N-1) = M(N, N-3)
    M(N+2, N)   = M(N, N-2)

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 3.0_dp)
      ! 使用 GD_basis_p4(4) 与 GD_basis_p4leftshift2(5)
      sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(5, xq + 2.0_dp)
    end do
    M(N+2, N+1) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1M1_1 + N1_N1M1_2 + N1_N1M1_3

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 3.0_dp)
      ! 使用 GD_basis_p4(4) 与 GD_basis_p4(4)
      sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(4, xq)
    end do
    M(N+2, N+2) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1_1 + N1_N1_2 + N1_N1_3

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 3.0_dp)
      ! 使用 GD_basis_p4(4) 与 GD_basis_p4rightshift2(3)
      sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(3, xq - 2.0_dp)
    end do
    M(N+2, N+3) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1P1_1 + N1_N1P1_2

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 3.0_dp)
      ! 使用 GD_basis_p4(4) 与 GD_basis_p4rightshift4(2)
      sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(2, xq - 4.0_dp)
    end do
    M(N+2, N+4) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1P2_1

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) + 3.0_dp)
      ! 使用 GD_basis_p4(4) 与 GD_basis_p4rightshift6(1)
      sum_val = sum_val + quad_weights(q) * phi_p4(4, xq) * phi_p4(1, xq - 6.0_dp)
    end do
    M(N+2, N+5) = (h/2.0_dp * sum_val)/2.0_dp

    ! 对于 M(N+3, *)（对应 φ10）
    M(N+3, N-1) = M(N, N-4)
    M(N+3, N)   = M(N, N-3)

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4leftshift4(5)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(5, xq + 4.0_dp)
    end do
    M(N+3, N+1) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1M2_1 + N1_N1M2_2

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4leftshift2(4)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(4, xq + 2.0_dp)
    end do
    M(N+3, N+2) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1M1_1 + N1_N1M1_2

    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4(3)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(3, xq)
    end do
    M(N+3, N+3) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1_1 + N1_N1_2
    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4rightshift2(2)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(2, xq - 2.0_dp)
    end do
    M(N+3, N+4) = (h/2.0_dp * sum_val)/2.0_dp + N1_N1P1_1
    sum_val = 0.0_dp
    do q = 1, numQuad
      xq = 0.5_dp * (quad_points(q) - 1.0_dp)
      ! 使用 GD_basis_p4(3) 与 GD_basis_p4rightshift4(1)
      sum_val = sum_val + quad_weights(q) * phi_p4(3, xq) * phi_p4(1, xq - 4.0_dp)
    end do
    M(N+3, N+5) = (h/2.0_dp * sum_val)/2.0_dp

  end subroutine assmblep4

  !===========================================================================
  ! 函数：phi_p4
  ! 功能：计算 P4 基函数值，index 取值 1~5，对应不同的多项式
  !===========================================================================
  real(dp) function phi_p4(index, x)
    implicit none
    integer, intent(in) :: index
    real(dp), intent(in) :: x
    select case(index)
    case(1)
      phi_p4 = (1.0_dp/384.0_dp) * (x+8.0_dp)*(x+6.0_dp)*(x+4.0_dp)*(x+2.0_dp)
    case(2)
      phi_p4 = -1.0_dp/96.0_dp * (x+6.0_dp)*(x+4.0_dp)*(x+2.0_dp)*(x-2.0_dp)
    case(3)
      phi_p4 = 1.0_dp/64.0_dp * (x+4.0_dp)*(x+2.0_dp)*(x-2.0_dp)*(x-4.0_dp)
    case(4)
      phi_p4 = -1.0_dp/96.0_dp * (x+2.0_dp)*(x-2.0_dp)*(x-4.0_dp)*(x-6.0_dp)
    case(5)
      phi_p4 = 1.0_dp/384.0_dp * (x-2.0_dp)*(x-4.0_dp)*(x-6.0_dp)*(x-8.0_dp)
    case default
      phi_p4 = 0.0_dp
    end select
  end function phi_p4

end program mass_matrix_assembly
