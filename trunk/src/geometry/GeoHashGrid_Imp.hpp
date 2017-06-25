// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        CGeoHashGrid<Real>::CGeoHashGrid() : m_nbCellsX(0), m_nbCellsY(0), m_nbCellsZ(0), m_nbCellsXY(0)
        {

            m_adOrig.resize(3);
            m_adDelta.resize(3);

        }

        template <typename Real>
        CGeoHashGrid<Real>::~CGeoHashGrid()
        {

        }

        template <typename Real>
        void CGeoHashGrid<Real>::reset()
        {

            CGeoContainer<CGeoCoordinate<Real>, Real>::reset();

            m_coordinateList.clear();
            m_coordinateListPtr.clear();

        }

        template <typename Real>
        void CGeoHashGrid<Real>::build()
        {

            try
            {

                // Get bounding box
                for (Integer i = 0; i < static_cast<Integer> (CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects.size()); ++i)
                    m_boundingBox.addCoordinate(CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i]);

                // Enlarge the box somehow - avoids division by zero along a dim.
                const Real C_EPS = std::numeric_limits<Real>::epsilon();

                // Compute dynamical size for bucket based on # points
                Real dCubicRoot = static_cast<Real> (pow(static_cast<Real> (CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects.size()), 0.333));

                m_nbCellsX = m_nbCellsY = m_nbCellsZ = static_cast<Integer>(dCubicRoot)+1;
                m_nbCellsXY = m_nbCellsX * m_nbCellsY;
                m_nbCells = m_nbCellsX * m_nbCellsY * m_nbCellsZ;

                m_coordinateListPtr.resize(m_nbCells + 1);

                m_adDelta = m_boundingBox.max() - m_boundingBox.min();

                if (m_adDelta[0] <= C_EPS)
                    m_nbCellsX = 1;

                if (m_adDelta[1] <= C_EPS)
                    m_nbCellsY = 1;

                if (m_adDelta[2] <= C_EPS)
                    m_nbCellsZ = 1;

                Real dDeltaMax = m_adDelta.maxCoeff();

                m_boundingBox.grow(C_EPS * dDeltaMax);

                // Find step size along each direction
                Real dRm1CellX = 1.0 / static_cast<Real>(m_nbCellsX);

                m_adDelta = (m_boundingBox.max() - m_boundingBox.min()) * dRm1CellX;

                m_adOrig = m_boundingBox.min();

                Real dRm1DeltaX = 1.0 / m_adDelta[0];
                Real dRm1DeltaY = 1.0 / m_adDelta[1];
                Real dRm1DeltaZ = 1.0 / m_adDelta[2];

                for (Integer ipass = 0; ipass < 2; ipass++)
                {

                    // Loop over pts - storing ahead

                    for (Integer i = 0; i < static_cast<Integer>(CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects.size()); ++i)
                    {

                        // Find in which bucket this point falls
                        Real x = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].x();
                        Real y = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].y();
                        Real z = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].z();

                        Integer imin = static_cast<Integer>((x - m_boundingBox.min().x()) * dRm1DeltaX);
                        Integer jmin = static_cast<Integer>((y - m_boundingBox.min().y()) * dRm1DeltaY);
                        Integer kmin = static_cast<Integer>((z - m_boundingBox.min().z()) * dRm1DeltaZ);

                        // Bound this value with correct values
                        imin = std::min(std::max(0, imin), static_cast<Integer> (m_nbCellsX - 1));
                        jmin = std::min(std::max(0, jmin), static_cast<Integer> (m_nbCellsY - 1));
                        kmin = std::min(std::max(0, kmin), static_cast<Integer> (m_nbCellsZ - 1));

                        // (imin,jmin,kmin) index - structured grid - cell indx
                        Integer indx = imin + jmin * m_nbCellsX + kmin * m_nbCellsXY;

                        if (ipass == 0)
                        {
                            // Increment count of pt falling into this cell - store ahead
                            m_coordinateListPtr[indx + 1]++;
                        }
                        else
                        {

                            // pointer to location where to store this pt
                            Integer istor = m_coordinateListPtr[indx];

                            // Store the pt
                            m_coordinateList[istor] = i;

                            // Incrt ptr to next available position
                            m_coordinateListPtr[indx]++;

                        }

                    }

                    if (ipass == 0)
                    {

                        // 1st PASS
                        // reshuffle  iptr_nodelist - Accumulate count
                        // iptr_nodelist[ie] : pts at the location where first entry 
                        // of pt should be stored 
                        for (Integer ie = 1; ie < m_nbCells + 1; ie++)
                        {
                            m_coordinateListPtr[ie] += m_coordinateListPtr[ie - 1];
                        }

                        // Total storage needed
                        Integer nstor = m_coordinateListPtr[m_nbCells];

                        // Allocate memory for nodelist
                        m_coordinateList.resize(nstor);

                    }
                    else
                    {

                        // 2nd PASS
                        // Finally reorder iptr_nodelist
                        for (Integer ie = m_nbCells; ie > 0; ie--)
                        {
                            m_coordinateListPtr[ie] = m_coordinateListPtr[ie - 1];
                        }

                        // Start pt
                        m_coordinateListPtr[0] = 0;

                    }

                }

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }

        }

        template <typename Real>
        void CGeoHashGrid<Real>::find(std::vector<Integer>& coordinateIds, CGeoCoordinate<Real>& aCoordinate, const Real aTolerance)
        {

            try
            {

                Real dRm1DeltaX = 1.0 / m_adDelta[0];
                Real dRm1DeltaY = 1.0 / m_adDelta[1];
                Real dRm1DeltaZ = 1.0 / m_adDelta[2];

                // Find cell where node falls - build a cube around the node
                Real dXMin = aCoordinate.x() - aTolerance;
                Real dYMin = aCoordinate.y() - aTolerance;
                Real dZMin = aCoordinate.z() - aTolerance;

                // Find in which bucket this point falls
                Integer imin = static_cast<Integer> ((dXMin - m_adOrig[0]) * dRm1DeltaX > std::numeric_limits<Integer>::min() ? (dXMin - m_adOrig[0]) * dRm1DeltaX : std::numeric_limits<Integer>::min());
                Integer jmin = static_cast<Integer> ((dYMin - m_adOrig[1]) * dRm1DeltaY > std::numeric_limits<Integer>::min() ? (dYMin - m_adOrig[1]) * dRm1DeltaY : std::numeric_limits<Integer>::min());
                Integer kmin = static_cast<Integer> ((dZMin - m_adOrig[2]) * dRm1DeltaZ > std::numeric_limits<Integer>::min() ? (dZMin - m_adOrig[2]) * dRm1DeltaZ : std::numeric_limits<Integer>::min());

                // Bound this value with correct values
                imin = std::min(std::max(0, imin), static_cast<Integer> (m_nbCellsX - 1));
                jmin = std::min(std::max(0, jmin), static_cast<Integer> (m_nbCellsY - 1));
                kmin = std::min(std::max(0, kmin), static_cast<Integer> (m_nbCellsZ - 1));

                Real dXMax = aCoordinate.x() + aTolerance;
                Real dYMax = aCoordinate.y() + aTolerance;
                Real dZMax = aCoordinate.z() + aTolerance;

                // Find in which bucket this point falls
                Integer imax = static_cast<Integer> ((dXMax - m_adOrig[0]) * dRm1DeltaX < std::numeric_limits<Integer>::max() ? (dXMax - m_adOrig[0]) * dRm1DeltaX : std::numeric_limits<Integer>::max());
                Integer jmax = static_cast<Integer> ((dYMax - m_adOrig[1]) * dRm1DeltaY < std::numeric_limits<Integer>::max() ? (dYMax - m_adOrig[1]) * dRm1DeltaY : std::numeric_limits<Integer>::max());
                Integer kmax = static_cast<Integer> ((dZMax - m_adOrig[2]) * dRm1DeltaZ < std::numeric_limits<Integer>::max() ? (dZMax - m_adOrig[2]) * dRm1DeltaZ : std::numeric_limits<Integer>::max());

                // Bound this value with correct values
                imax = std::min(std::max(0, imax), static_cast<Integer> (m_nbCellsX - 1));
                jmax = std::min(std::max(0, jmax), static_cast<Integer> (m_nbCellsY - 1));
                kmax = std::min(std::max(0, kmax), static_cast<Integer> (m_nbCellsZ - 1));

                // Loop over cells Intersecting node Bounding Box
                for (int i = imin; i <= imax; ++i)
                {

                    for (int j = jmin; j <= jmax; ++j)
                    {

                        Integer j_off = j * m_nbCellsX;

                        for (int k = kmin; k <= kmax; ++k)
                        {

                            Integer ie = i + j_off + k * m_nbCellsXY;

                            // Ptrs to start and end nodes in this cell
                            Integer ip_start = m_coordinateListPtr[ie];
                            Integer ip_end = m_coordinateListPtr[ie + 1];

                            // Loop over points in these cells
                            for (Integer iptr = ip_start; iptr < ip_end; iptr++)
                            {

                                Integer aCoordinateIndex = m_coordinateList[iptr];

                                Real dist = (CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[aCoordinateIndex] - aCoordinate).norm();

                                if (dist <= aTolerance)
                                    coordinateIds.push_back(CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjectIds[aCoordinateIndex]);

                            }

                        }

                    }

                }
            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }

        }

    }

}
