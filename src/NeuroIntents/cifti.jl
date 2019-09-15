# DISCLAIMER!!! Everything in this file should be considered pre-alpha and subject to change
"""
    DenseConnectivity

AppliesToMatrixDimension 0: brain models
AppliesToMatrixDimension 1: brain models

This file type represents connectivity between points in the brain. A row is
the connectivity from a single vertex or voxel in the mapping that applies
along the second dimension, to all vertices and voxels along the first
dimension. This specification of “from” and “to” is not intended to imply that
the data is always directed, but to establish a convention so that directed
data is stored consistently, and to ensure that interactive usage loads rows
from the matrix, in order to maximize responsiveness. Note that this type
can have a single mapping apply to both dimensions, but can also have separate,
different mappings for rows and columns, for instance only containing
connectivity from left cortex to cerebellum.
"""
struct DenseConnectivity end
const ConnDense = DenseConnectivity()
intentcode(::DenseConnectivity) = 3001
intentname(::DenseConnectivity) = "ConnDense"
intent2ext(::DenseConnectivity) = "dconn.nii"

"""
    ConnectivityDenseDataSeries

AppliesToMatrixDimension 0: series
AppliesToMatrixDimension 1: brain models

This file type represents data points in a series for every vertex and voxel in
the mapping. A row is a complete data series, for a single vertex or voxel in
the mapping that applies along the second dimension. A dataseries is often a
timeseries, but it can also represent other data types such as a series of
sampling depths along the surface normal from the white to pial surface. It
retains the ‘t’ in dtseries from CIFTI-1 naming conventions.
"""
struct ConnectivityDenseDataSeries end
const ConnDenseSeries = ConnectivityDenseDataSeries()
intentcode(::ConnectivityDenseDataSeries) = 3002
intentname(::ConnectivityDenseDataSeries) = "ConnDenseSeries"
intent2ext(::ConnectivityDenseDataSeries) = ".dtseries.nii"

"""
    ParcellatedConnectivity

AppliesToMatrixDimension 0: parcels
AppliesToMatrixDimension 1: parcels

This file type represents connectivity between areas or parcels of the brain.
Similarly to Dense Connectivity, a row is the connectivity from a single parcel
in the mapping that applies along the 14 second dimension, to all parcels along
the first dimension. Note that this type can have a single mapping apply to
both dimensions, but can also be asymmetric, for example if a parcellated
connection was from a subset of cortical areas to all of them (e.g. injection
sites for tracer data).
"""
struct ParcellatedConnectivity end
const ConnParcels = ParcellatedConnectivity()
intentcode(::ParcellatedConnectivity) = 3003
intentname(::ParcellatedConnectivity) = "ConnParcels"
intent2ext(::ParcellatedConnectivity) = ".pconn.nii"

"""
    ParcellatedDataSeries

AppliesToMatrixDimension 0: series
AppliesToMatrixDimension 1: parcels

This file type represents data points in a series for areas of the brain.
Similarly to Dense Data Series, a row is a complete data series (e.g. ICA
component timeseries), for a single parcel in the mapping that applies along
the second dimension. This is called ptseries by analogy with dtseries.
"""
struct ParcellatedDataSeries end
const ConnParcelSries = ParcellatedDataSeries()
intentcode(::ParcellatedDataSeries) = 3004
intentname(::ParcellatedDataSeries) = "ConnParcelSries"
intent2ext(::ParcellatedDataSeries) = ".ptseries.nii"

"""
    ConnectivityDenseScalars

File extension: .dscalar.nii, others - see “Specializations of Scalar Maps”
AppliesToMatrixDimension 0: scalars
AppliesToMatrixDimension 1: brain models
This file type stores named scalar maps across vertices and voxels, like a GIFTI “functional” file
(aka a metric file) that can accommodate multiple surfaces, and voxels. It is often used to store
task fMRI activity maps in the standard grayordinates space or myelin maps from both cortical
hemispheres. A row is a single value from each named map, for one vertex or voxe
"""
struct ConnectivityDenseScalars end
const ConnDenseScalar = ConnectivityDenseScalars()
intentcode(::ConnectivityDenseScalars) = 3006
intentname(::ConnectivityDenseScalars) = "ConnDenseScalar"
FileIO.file_extension(::ConnectivityDenseScalars) = ".dscalar.nii"

"""
    ConnectivityDenseLabel

AppliesToMatrixDimension 0: labels
AppliesToMatrixDimension 1: brain models
This file type stores named label maps across vertices and voxels, similarly to
Dense Scalar, but with label keys in the data matrix instead of
continuous-valued data. Each label map has its own separate label table defined
within the CIFTI XML. A row is a single label value from each named label map,
for one vertex or voxel.
"""
struct ConnectivityDenseLabel end
const ConnDenseLabel = ConnectivityDenseLabel()
intentcode(::ConnectivityDenseLabel) = 3007
intentname(::ConnectivityDenseLabel) = "ConnDenseLabel"
intent2ext(::ConnectivityDenseLabel) = ".dlabel.nii"


"""
    ConnectivityParcellatedScalar

AppliesToMatrixDimension 0: scalars

This file type stores named scalar maps across parcels. A row is a single
scalar from each map, for one parcel. For example, it could store task activity
measures in parcels.
"""
struct ConnectivityParcellatedScalar end
const ConnParcelScalr = ConnectivityParcellatedScalar()
intentcode(::ConnectivityParcellatedScalar) = 3008
intentname(::ConnectivityParcellatedScalar) = "ConnParcelScalr"
intent2ext(::ConnectivityParcellatedScalar) = ".pscalar.nii"


"""
    ParcellatedDenseConnectivity

AppliesToMatrixDimension 0: brain models
AppliesToMatrixDimension 1: parcels

This file type stores connectivity from parcels to vertices and/or voxels. A
row is the connectivity from one parcel to all vertices and voxels.
"""
struct ParcellatedDenseConnectivity end
const ConnParcelDense = ParcellatedDenseConnectivity()
intentcode(::ParcellatedDenseConnectivity) = 3009
intentname(::ParcellatedDenseConnectivity) = "ConnParcelDense"
intent2ext(::ParcellatedDenseConnectivity) = ".pdconn.nii"

"""
    DenseParcellatedConnectivity

AppliesToMatrixDimension 0: parcels
AppliesToMatrixDimension 1: brain models

This file type stores connectivity from vertices and voxels to parcels. A row
is the connectivity from one vertex or voxel to all parcels.
"""
struct DenseParcellatedConnectivity end
const ConnDenseParcel = DenseParcellatedConnectivity()
intentcode(::DenseParcellatedConnectivity) = 3010
intentname(::DenseParcellatedConnectivity) = "ConnDenseParcel"
intent2ext(::DenseParcellatedConnectivity) = ".dpconn.nii"

"""
    ParcellatedConnectivityScalar

This file type stores connectivity between parcels under named conditions (e.g.
frequency bands, different subjects). A row is the connectivity from one parcel
in the mapping that applies along the second dimension to all parcels in the
first dimension at one named condition (e.g. frequency bands) in the third
dimension.
"""
struct ParcellatedConnectivityScalar end
const ConnPPSc = ParcellatedConnectivityScalar()
intentcode(::ParcellatedConnectivityScalar) = 3013
intentname(::ParcellatedConnectivityScalar) = "ConnPPSc"
intent2ext(::ParcellatedConnectivityScalar) = "pconnscalar.nii"

# TODO
"""
    DenseFiberFans

AppliesToMatrixDimension 0: scalars
AppliesToMatrixDimension 1: brain models

Map Layout:
For brainordinate index B:
Index (0, B): X coordinate
Index (1, B): Y coordinate
Index (2, B): Z coordinate

For fiber F, starting at 0:
Index (7 * F + 3, B): Mean of fiber strength
Index (7 * F + 4, B): Standard deviation of fiber strength
Index (7 * F + 5, B): Theta angle of fiber bingham
Index (7 * F + 6, B): Phi angle of fiber bingham
Index (7 * F + 7, B): Dispersion parameter ka of fiber bingham
Index (7 * F + 8, B): Dispersion parameter kb of fiber bingham
Index (7 * F + 9, B): Psi angle of fiber bingham

This file is used to store fiber fan bingham models (either of fiber uncertainty—e.g. from
bedpost—or fiber fanning), for the purpose of displaying probabilistic tractography data
(Sotiropoulos et al., 2012)
"""
struct DenseFiberFan end
const FiberFan = DenseFiberFan()
# note that Dense fiber fan is technically not a unique intent
intentcode(::DenseFiberFan) = 3002
intentname(::DenseFiberFan) = "ConnDenseSeries"
intent2ext(::DenseFiberFan) = "dfan.nii"

"""
    DenseFanSamples

AppliesToMatrixDimension 0: scalars
AppliesToMatrixDimension 1: scalars
AppliesToMatrixDimension 2: brain models
Map Layout:
For brainordinate index B:
Index (0, 0, B): X coordinate
Index (1, 0, B): Y coordinate
Index (2, 0, B): Z coordinate
Other indexes in (x, 0, B) are unused
For sample S, starting at 0:
For fiber F, starting at 0:
Index (6 * F, S + 1, B): Fan strength
Index (6 * F + 1, S + 1, B): Theta angle of fan bingham
Index (6 * F + 2, S + 1, B): Phi angle of fan
Index (6 * F + 3, S + 1, B): Dispersion parameter ka of fan bingham
Index (6 * F + 4, S + 1, B): Dispersion parameter kb of fan bingham
Index (6 * F + 5, S + 1, B): Psi angle of fan bingham

This file is used for a similar purpose as dense fiber samples, but is used
when fiber fanning is modeled directly.
"""
struct DenseFanSamples end
const FanSamples = DenseFanSamples()
intentcode(::DenseFanSamples) = 3000
intentname(::DenseFanSamples) = "ConnUnkown"
intent2ext(::DenseFanSamples) = ".dfansamp.nii"


struct ConnectivityUnkown end
const ConnUnkown = ConnectivityUnkown()
intentcode(::ConnectivityUnkown) = 3000
intentname(::ConnectivityUnkown) = "ConnUnkown"

function ciftiintent(i::Integer)
    if i == 3000
        return ConnUnkown
    elseif i == 3001
        return ConnDense
    elseif i == 3002
        return ConnDenseSeries
    elseif i == 3003
        return ConnParcels
    elseif i == 3004
        return ConnParcelSeries
    elseif i == 3006
        return ConnDenseScalar
    elseif i == 3007
        return ConnDenseLabel
    elseif i == 3008
        return ConnParcelScalr
    elseif i == 3009
        return ConnParcelDense
    elseif i == 3010
        return ConnDenseParcel
    elseif i == 3013
        return ConnPPSc
    end
end

function ext2intent(ext::Vector{String})
    if length(ext) > 1
        # note the extra file extension types specify that they shouldn't be
        # gzipped so this should be as far as we need to parse the extensions
        # on gz files.
        if last(ext) == "gz"
            return GZip
        elseif ext[end-1] == "dconn"
            return ConnDense
        elseif ext[end-1] == "dtseries"
            return ConnDenseSeries
        elseif ext[end-1] == "pconn"
            return ConnParcels
        elseif ext[end-1] == "ptseries"
            return ConnParcelSries
        elseif ext[end-1] == "dscalar"
            return ConnDenseScalar
        elseif ext[end-1] == "dlabel"
            return ConnDenseLabel
        elseif ext[end-1] == "pscalar"
            return ConnParcelScalr
        elseif ext[end-1] == "pdconn"
            return ConnParcelDense
        elseif ext[end-1] == "dpconn"
            return ConnDenseParcel
        elseif ext[end-1] == "pconnscalar"
            return ConnPPSc
        elseif ext[end-1] == "dfan"
            return FiberFan
        elseif ext[end-1] == "dfansamp"
            return FanSamples
        end
    else
    end
end

