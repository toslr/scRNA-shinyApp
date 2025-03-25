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

// Initialize collapsible sidebar sections when the document is ready
$(document).ready(function() {
  // Set up click handlers for collapsible sections
  $('.collapsible-section .section-header').on('click', function() {
    var section = $(this).parent();
    section.toggleClass('collapsed');
    
    // Store the state in localStorage to remember it between sessions
    var sectionId = section.index();
    localStorage.setItem('section_' + sectionId + '_collapsed', section.hasClass('collapsed'));
  });
  
  // Initialize sections based on saved state (if any)
  $('.collapsible-section').each(function(index) {
    var defaultCollapsed = [2,3,4,5,6];
    var isCollapsed = defaultCollapsed.includes(index);
    if (isCollapsed) {
      $(this).addClass('collapsed');
    }
  });
  
  // Existing scroll function from your original script.js
  window.scrollToSection = function(sectionId) {
    var section = document.getElementById(sectionId);
    if (section) {
      section.scrollIntoView({ behavior: 'smooth' });
    }
  };
});