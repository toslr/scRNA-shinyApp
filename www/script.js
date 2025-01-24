function scrollToSection(sectionId) {
  const element = document.getElementById(sectionId);
  const topBar = document.getElementById('topbar');
  const topBarHeight = topBar.offsetHeight;
  const elementPosition = element.getBoundingClientRect().top + window.pageYOffset;
  
  // Add 20px extra padding
  const offsetPosition = elementPosition - topBarHeight - 20;

  window.scrollTo({
    top: offsetPosition,
    behavior: 'smooth'
  });
}