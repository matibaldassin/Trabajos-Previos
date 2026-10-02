program vtk_local_polarization_strain
implicit none

integer, parameter :: celdas = 21952
integer, parameter :: bsite = celdas

real(8) :: px, py, pz, p, identidad
real(8) :: strain31, strain32, strain33
real(8) :: xx, yy, zz
integer :: i, j, k, n, centros
integer :: levcfg, imcon, natoms, lines_after_coordinates
integer :: point_index
character(20) :: atom
character(200) :: header
character(20) :: species
character(500) :: csv_line

open(1, file="CONFIG", status="old", action="read")
open(2, file="entrada.dat", status="old", action="read")
open(10, file="salida-no-angles-deltaP.vtk", status="replace", action="write")
open(22, file="ID", status="replace", action="write")

write(10,'(A)') "# vtk DataFile Version 2.0"
write(10,'(A)') "polydata legacy vtk file with point scalar data"
write(10,'(A)') "ASCII"
write(10,'(A)') "DATASET POLYDATA"
write(10,'(A, I6, A)') "POINTS", bsite, "  double"

read(1,*)
read(1,*) levcfg, imcon, natoms

do i = 1, 3
   read(1,*)
enddo

if (levcfg.eq.0) then
   lines_after_coordinates = 0
else if (levcfg.eq.1) then
   lines_after_coordinates = 1
else if (levcfg.eq.2) then
   lines_after_coordinates = 2
else
   stop "Unsupported CONFIG levcfg value"
endif

do i = 1, celdas
   read(1,*) atom, n
   read(1,*) xx, yy, zz
   do j = 1, lines_after_coordinates
      read(1,*)
   enddo
   
   read(1,*) atom, n
   read(1,*) xx, yy, zz
   do j = 1, lines_after_coordinates
      read(1,*)
   enddo

   write(10,'(3(F12.6,3x))') xx, yy, zz

   if (trim(atom).eq."Zr") then
      write(22,*) 1
   else if (trim(atom).eq."Ti") then
      write(22,*) 2
   else
      write(22,*) 3
   endif

   do k = 3, 10
      read(1,*)
      read(1,*)
      do j = 1, lines_after_coordinates
         read(1,*)
      enddo
   enddo
enddo

close(1)
close(22)

write(10,'(A, I6, 1x, I6)') "VERTICES", bsite, 2*bsite
do i = 1, bsite
   point_index = i - 1
   write(10,'(I1, 1x, I6)') 1, point_index
enddo

write(10,'(A, I6)') "POINT_DATA", bsite

write(10,'(A)') "SCALARS polar_module double"
write(10,'(A)') "LOOKUP_TABLE default"
do i = 1, bsite
   read(2,*) n, px, py, pz, p
   write(10,'(F12.6)') p
enddo
close(2)

open(22, file="ID", status="old", action="read")
write(10,'(A)') "SCALARS identidad double"
write(10,'(A)') "LOOKUP_TABLE default"
do i = 1, bsite
   read(22,*) identidad
   write(10,'(F3.1)') identidad
enddo
close(22)

open(2, file="entrada.dat", status="old", action="read")
write(10,'(A)') "VECTORS polar_vector double"
do i = 1, bsite
   read(2,*) n, px, py, pz, p
   write(10,'(3(F12.6,2x))') px, py, pz
enddo
close(2)
close(10)

end program vtk_local_polarization_strain

