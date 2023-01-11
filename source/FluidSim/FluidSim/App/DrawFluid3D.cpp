#include "DrawFluid3D.h"
#include "glad/glad.h"

DrawFluid3D::DrawFluid3D(const Fluid3D& inFluid)
{
	FromFluid(inFluid);
}

void DrawFluid3D::FromFluid(const Fluid3D& inFluid)
{
	mParticlePoints.clear();
	mParticlePoints.reserve(inFluid.GetNumParticles());

	std::vector<Particle3D>::const_iterator pItr;

	for (pItr = inFluid.GetParticles().begin(); pItr < inFluid.GetParticles().end(); ++pItr)
	{
		mParticlePoints.push_back(pItr->GetPosition());
	}
}

void DrawFluid3D::UpdateOpenGLBuffers()
{

}

